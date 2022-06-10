
#%% Import Libraries
import numpy as np
import pybullet as p
import pybullet_data
import pinocchio as pin
import math

import crocoddyl
import time
import cv2 as cv
from projection_modified import pointProject, projectMatrix
import matplotlib.pyplot as plt
pin.switchToNumpyArray()

#%% functions
class ActionModelRobot(crocoddyl.ActionModelAbstract):
    def __init__(self, state, actuationModel, costModel):
        crocoddyl.ActionModelAbstract.__init__(self, state, state.nv )
        self.actuation = actuationModel
        self.costs = costModel
        
    def init_robot_sys(self,robot_sys, nr = 1):
        self.robot_sys = robot_sys
        self.Du = robot_sys.Du
        self.Dx = robot_sys.Dx
        self.Dr = nr
        
    def calc(self, data, x, u):
        q, v = x[:self.state.nq], x[-self.state.nv:]
        pin.forwardKinematics(self.state.pinocchio, data.pinocchio, q, v)
        pin.updateFramePlacements(self.state.pinocchio, data.pinocchio)
        pin.computeAllTerms(self.state.pinocchio, data.pinocchio, q, v)
        self.costs.calc(data.costs,x,u)
        data.cost=data.costs.cost
        data.xnext = self.robot_sys.step(x,u)
        
    def calcDiff(self, data, x, u, recalc = True):
        if recalc:
            self.calc(data, x, u)
        nq, nv = self.state.nq, self.state.nv
        q, v = x[:nq], x[-nv:]
        pin.computeAllTerms(self.state.pinocchio, data.pinocchio, q, v)
        pin.computeRNEADerivatives(self.state.pinocchio, data.pinocchio, q, v, pin.utils.zero(nv)) # Zero for acceleration because of signle integrator system
        self.costs.calcDiff(data.costs, x, u)
        A, B = self.robot_sys.compute_matrices(x,u)
        data.Fx = A.copy()
        data.Fu = B.copy()
        
    def createData(self):
        data = ActionDataRobot(self)
        return data

class ActionDataRobot(crocoddyl.ActionDataAbstract):
    def __init__(self, model):
        crocoddyl.ActionDataAbstract.__init__(self,model)
        self.pinocchio = pin.Model.createData(model.state.pinocchio)
        self.multibody = crocoddyl.DataCollectorMultibody(self.pinocchio)
        self.actuation = model.actuation.createData()
        self.costs = model.costs.createData(self.multibody)
        self.costs.shareMemory(self)
        self.Minv = None

class Single_integ():
    def __init__(self, rmodel, dt = 0.01):
        self.dt = dt
        self.rmodel = rmodel
        self.data = rmodel.createData()
        self.Dx = rmodel.nq + rmodel.nv 
        self.Du = rmodel.nv
           
    def compute_matrices(self,x,u):
        A = np.zeros((self.Dx,self.Dx))
        B = np.zeros((self.Dx,self.Du))
        A[:int(self.Dx/2),:int(self.Dx/2)] = np.eye(int(self.Dx/2))
        A[:int(self.Dx/2),int(self.Dx/2):] = np.eye(int(self.Dx/2)) * self.dt
        B[int(self.Dx/2):,:] = np.eye(self.Du)
        self.A, self.B = A, B
        return A,B
        
    def step(self, x, u):
        self.compute_matrices(x,u)
        x_next = self.A @ x + self.B @ u
        return x_next
    
    def rollout(self,x0, us):
        x_cur = x0.copy()
        xs = [x_cur]
        T = len(us)
        for i in range(T):
            x_cur = self.step(x_cur, us[i])
            xs += [x_cur]
        return np.array(xs)

class ObstacleAvoidance(crocoddyl.CostModelAbstract):
    def __init__(self, state, frameId, Q2, dist, ThetaX, ThetaY, target, traj, traj_weights, activation=None, nu=None):
        activation = activation if activation is not None else crocoddyl.ActivationModelQuad(3)
        if nu is None:
            crocoddyl.CostModelAbstract.__init__(self, state, activation)
        else:
            crocoddyl.CostModelAbstract.__init__(self, state, activation, nu)
        # self.obsPos = obsPos
        self.frameId = frameId
        self.dist = dist
        self.pitch = ThetaX
        self.yaw = ThetaY
        self.target = target
        self.Q2 = -Q2
        self.traj = traj
        self.traj_weights = traj_weights
    def calc(self, data, x, u): 
        nq, nv = self.state.nq, self.state.nv
        q, v = x[:nq], x[-nv:]
        pin.computeAllTerms(self.state.pinocchio, data.shared.pinocchio, q, v)
        pin.updateFramePlacements(self.state.pinocchio, data.shared.pinocchio)
        ee_pos = data.shared.pinocchio.oMf[self.frameId].translation
        c = 0
        for repPoint, w in zip(self.traj, self.traj_weights):
            obs_proj, _ = pointProject(repPoint, self.yaw, self.pitch, self.dist, self.target)
            ee_proj, _  = pointProject(ee_pos, self.yaw, self.pitch, self.dist, self.target)
            d = ee_proj - obs_proj
            c += w*d.transpose()@ self.Q2 @ d
        # data.r = data.shared.pinocchio.oMf[self.frameId].translation - self.obsPos
        # self.activation.calc(data.activation, data.r)
        data.cost =  c

    def calcDiff(self, data, x, u): 
        nq, nv = self.state.nq, self.state.nv
        q, v = x[:nq], x[-nv:]
        pin.computeAllTerms(self.state.pinocchio, data.shared.pinocchio, q, v)
        pin.updateFramePlacements(self.state.pinocchio, data.shared.pinocchio)
        self.activation.calcDiff(data.activation, data.r)
        J = pin.computeFrameJacobian(self.state.pinocchio, data.shared.pinocchio, q, self.frameId, pin.LOCAL_WORLD_ALIGNED)
        J_trans = J[:3,:]
        dq = np.zeros([nq,])
        dqq = np.zeros([nq,nq])
        ee_pos = data.shared.pinocchio.oMf[self.frameId].translation
        for repPoint, w in zip(self.traj, self.traj_weights):
            obs_proj, R = pointProject(repPoint, self.yaw, self.pitch, self.dist, self.target)
            ee_proj, _  = pointProject(ee_pos, self.yaw, self.pitch, self.dist, self.target)
            d = ee_proj - obs_proj
            dq += 2 * w*J_trans.transpose() @ R.transpose() @ self.Q2 @ d
            dqq += 2 * w*J_trans.transpose() @ R.transpose() @ self.Q2 @ R @ J_trans
        # data.Rx = np.hstack([J[:3,:], pin.utils.zero((self.activation.nr, self.state.nv))])
        data.Lx[:nq] = dq
        data.Lxx [:nq,:nq] = dqq
        # data.Lxx = np.dot(data.Rx.T, np.dot(data.activation.Arr, data.Rx))

class ObstacleAvoidance2(crocoddyl.CostModelAbstract):
    def __init__(self, state, frameId,  eigen_vec, dist, ThetaX, ThetaY, target, Q2, dmax, traj_weights = None, traj=None, nu=None):
        
        activation = crocoddyl.ActivationModelQuad(0)

        if nu is None:
            crocoddyl.CostModelAbstract.__init__(self, state, activation)
        else:
            crocoddyl.CostModelAbstract.__init__(self, state, activation, nu)
        self.traj = traj
        self.frameId = frameId
        # self.Q2 = Q2
        if traj_weights is None:
            self.traj_weights = np.ones(len(self.traj))
        else:
            self.traj_weights = traj_weights

        self.eigen_vec = eigen_vec
        self.dist = dist
        self.pitch = ThetaX
        self.yaw = ThetaY
        self.target = target
        self.Q2 = Q2
        self.dmax = dmax

    def calc(self, data, x, u): 
        nq, nv = self.state.nq, self.state.nv
        q, v = x[:nq], x[-nv:]
        pin.computeAllTerms(self.state.pinocchio, data.shared.pinocchio, q, v)
        pin.updateFramePlacements(self.state.pinocchio, data.shared.pinocchio)
        xd = data.shared.pinocchio.oMf[self.frameId].translation.transpose()
        c = 0

        for repPoint, w,  in zip(self.traj, self.traj_weights):
            d1,_ = pointProject(xd, self.yaw, self.pitch, self.dist, self.target)
            d2,_ = pointProject(repPoint, self.yaw, self.pitch, self.dist, self.target)
            d = d1 - d2
            # d, _ = pointProject(d, self.ThetaX, self.ThetaY, self.ThetaZ, self.camPos)
            # d = self.eigen_vec.transpose() @ d
            # d = w*d
            d_l1 = np.absolute(d)
            f1 = 0
            f2 = 0
            f3 = 0
            if d_l1[0] < self.dmax:
                f1 = self.dmax - d_l1[0]
            if d_l1[1] < self.dmax:
                f2 = self.dmax - d_l1[1]
            f = np.array([[f1], [f2]])

            c = c + w* f.transpose() @ self.Q2 @ f  
        data.cost = c[0,0]
        # test = 1

    def calcDiff(self, data, x, u): 
        nq, nv = self.state.nq, self.state.nv
        q, v = x[:nq], x[-nv:]
        pin.computeAllTerms(self.state.pinocchio, data.shared.pinocchio, q, v)
        pin.updateFramePlacements(self.state.pinocchio, data.shared.pinocchio)
        #self.activation.calcDiff(data.activation, data.r)
        J = pin.computeFrameJacobian(self.state.pinocchio, data.shared.pinocchio, q, self.frameId, pin.LOCAL_WORLD_ALIGNED)
        xd = data.shared.pinocchio.oMf[self.frameId].translation.transpose()
        Jx = np.zeros((2,1))
        Hxx = np.zeros((2,2))
        R = projectMatrix(self.yaw,self.pitch)
        
        R = R[[0,1],:]

        for repPoint, w in zip(self.traj, self.traj_weights):
            d1,_ = pointProject(xd, self.yaw, self.pitch, self.dist, self.target)
            d2,_ = pointProject(repPoint, self.yaw, self.pitch, self.dist, self.target)
            d = d1 - d2
            # d, _ = pointProject(d, self.ThetaX, self.ThetaY, self.ThetaZ, self.camPos)
            # d =  d
            # d = d
            # d = (xd-repPoint)
            # d, _ = pointProject(d, self.ThetaX, self.ThetaY, self.ThetaZ, self.camPos)
            # d = self.eigen_vec.transpose() @ d
            # d = w*d 
            d_l1 = np.absolute(d)
            j1 = 0
            j2 = 0
            f1 = 0
            f2 = 0
        
            if d_l1[0] < self.dmax:
                f1 = self.dmax - d_l1[0]
                if repPoint[0] < xd[0]: j1 = -1
                else: j1 = 1
            
            if d_l1[1] < self.dmax:
                f2 = self.dmax - d_l1[1]
                if repPoint[1] < xd[1]: j2 = -1
                else: j2 = 1
            
            f = np.array([[f1], [f2]])
            Jac = np.diag([j1, j2])

            Jx += 2*Jac.transpose() @ self.Q2 @ f * w
            Hxx += Jac.transpose() @ self.Q2 @ Jac * w
    

        test = 1             
        #Lq = J[0:3,:].transpose() @ R[0,:].reshape([3,1]) @ Jx[0] + J[0:3,:].transpose() @ R[1,:].reshape([3,1]) @ Jx[1] 
        
        data.Lx = np.vstack([J[:3,:].transpose() @ R.transpose() @ Jx, np.zeros((nv,1))]) 
        Lxx = np.zeros((nq+nv, nq+nv))
        Lxx[0:nq,0:nq] = J[:3,:].transpose() @  R.transpose() @ Hxx @ R @ J[:3,:]

        data.Lxx = Lxx
        """
        ddq_dq = np.zeros([nv,nv],dtype=np.float64)
        J1 = J.copy()
        for i in range (0,nv):
            q2 = q.copy()
            dx = 1e-3
            q2[i] =  q2[i] + dx
            pin.computeAllTerms(self.state.pinocchio, data.shared.pinocchio, q2, v)
            pin.computeRNEADerivatives(self.state.pinocchio, data.shared.pinocchio, q2, v, v)
            J2 = pin.getFrameJacobian(self.state.pinocchio,data.shared.pinocchio,self.frameId,pin.LOCAL_WORLD_ALIGNED)
            dJ_dq_i = ((J2-J1)/dx)
            ddq_dq[:,i]  =dJ_dq_i[0:3,:].transpose() @ R[0,:].reshape([3,1]) @ Jx[0] + dJ_dq_i[0:3,:].transpose() @ R[1,:].reshape([3,1]) @ Jx[1] 
        Lxx = np.zeros((nq+nv, nq+nv))
        Lxx[0:nq,0:nq] = ddq_dq
        data.Lxx = Lxx
        """
        
def get_joint_limits(robot_id, joints):
    lower_limits = []
    upper_limits = []
    for i in range(len(joints)):
        info = p.getJointInfo(robot_id, joints[i])
        lower_limits += [info[8]]
        upper_limits += [info[9]]
    limits = np.vstack([lower_limits, upper_limits])
    return limits

def check_limits(q, joint_limits):
    for i in range(q.shape[0]):
        if q[i] < joint_limits[0,i] - 1e-10 or q[i] > joint_limits[1,i] + 1e-10:
            raise ValueError("Joint limits are violated at joint {}. q = {}".format(i,q))
            
#%% parameters
URDF = '/home/student/simulation_example/URDF/frankaemika_new/panda_arm.urdf'
ee = 'panda_grasptarget'
hand_idx = 7 # I think for pybullet
dt = 0.1
duration = 7 #total time that you expect the robot finish the task
T = int(duration/dt)
pos_goal = np.array([0.5,0.1,0.1])
obsPos = np.array([0.5,-0.1,0.1])
max_force = 500
control_mode = p.POSITION_CONTROL
movie_flag = 0
pixelWidth = 1280
pixelHeight = 720
camTargetPos = [0, 0, 0]
camDistance = 2
camPos = [camDistance,0,0]
pitch = -20
roll = 0
yaw = 90
upAxisIndex = 2 
movie_save = "./movies/Panda_Legible(pitch = {}, roll = {}, yaw = {}, TargetPos = {}, obstaclePos = {}).avi".format(pitch,roll,yaw, camTargetPos, obsPos)

#%% Initialize the robot model
rmodel = pin.buildModelFromUrdf(URDF)
rdata = rmodel.createData()
eeIdx = rmodel.getFrameId(ee)
state = crocoddyl.StateMultibody(rmodel)
actuation = crocoddyl.ActuationModelFull(state)
q0 = np.array([0, -np.pi/3, 0, -3*np.pi/4, 0, 3*np.pi/4, np.pi/4, 0.02, 0.02])
x0 = np.concatenate([q0, np.zeros([state.nv,])])
dmax = 0.3
#%% Building obstacle trajectories

# start by updating position of the frames
pin.forwardKinematics(rmodel, rdata, q0, np.zeros([state.nv]))
pin.updateFramePlacements(rmodel, rdata)
fwd_x0T = rdata.oMf[eeIdx].translation

# build trajectory of repulsive points
M = 30              # number of repulsive seed points
traj = np.zeros((M, 3))
traj[:,0] = np.linspace(fwd_x0T[0], obsPos[0], num=M)
traj[:,1] = np.linspace(fwd_x0T[1], obsPos[1], num=M)
traj[:,2] = np.linspace(fwd_x0T[2], obsPos[2], num=M)
obs_poses = traj
#%% look for eigenvalues vectors to build

pos_goal_proj, R = pointProject(pos_goal.transpose(), yaw, pitch, camDistance, camTargetPos)
obs_proj,_ = pointProject(obsPos.transpose(), yaw, pitch, camDistance, camTargetPos)
R = R[[0,1],:]
X = np.array([pos_goal_proj, obs_proj]).transpose()
X = X-np.mean(X, axis=1)
CovMat = X @ X.transpose()
eigVal, eigVec = np.linalg.eig(CovMat)

# sorts the eigenvecors in descending order with respect to eigenvalues
sortIdx = np.argsort(eigVal)
sortIdx = sortIdx[::-1]
eigVec = eigVec[:,sortIdx]

#%% Creating the weight matrix for the repuslive cost
Q2 = np.diag([5, 1])
Q2 = eigVec @ Q2 @ eigVec.transpose()
traj_weights = np.ones((len(traj),2))*1
traj_weights[:,0] = np.linspace(2,0,len(traj))**2
traj_weights[:,1] = np.linspace(0.2,0,len(traj))**2

traj_weights =np.ones(len(traj))*2
#traj_weights =np.linspace(4, 0, len(traj))**2

#%% Cost function
# boundries
xlb = np.concatenate((rmodel.lowerPositionLimit, -100 * np.ones((state.nv)))) #velocity boundry is not from the real robot
xub = np.concatenate((rmodel.upperPositionLimit, 100 *  np.ones((state.nv)))) #velocity boundry is not from the real robot
bounds = crocoddyl.ActivationBounds(xlb, xub, 1.)
xLimitResidual = crocoddyl.ResidualModelState(state, pin.utils.zero([rmodel.nq + rmodel.nv]), actuation.nu)
xLimitActivation = crocoddyl.ActivationModelQuadraticBarrier(bounds)
limitCost = crocoddyl.CostModelResidual(state, xLimitActivation, xLimitResidual)

# goal tracking
goalTrackingres = crocoddyl.ResidualModelFrameTranslation(state, eeIdx, pos_goal)
goalTrackingActivation = crocoddyl.ActivationModelWeightedQuad(np.matrix(np.array([1] * 3)).T)
goalTrackingCost = crocoddyl.CostModelResidual(state, goalTrackingActivation, goalTrackingres)

#ObstacleAvoidance 
# obsAvoidres = crocoddyl.ResidualModelFrameTranslation(state, eeIdx, obsPos)
# obsAvoidActivation = crocoddyl.ActivationModelWeightedQuad(-Q2)
#obsAvoidCost = ObstacleAvoidance(state, eeIdx, activation = obsAvoidActivation, obsPos=obsPos)
obsAvoidCost = ObstacleAvoidance(state, eeIdx, Q2, dist = camDistance, ThetaX=pitch, ThetaY=yaw, target = camTargetPos, traj=traj, traj_weights=traj_weights)
# default state ||x-x_def||^2 (x regulator)
defaultState = x0
stateWeights = np.array([.1] * state.nq + [100] * state.nv)
stateWeightsTerm = np.array([.1] * state.nq + [100] * state.nv)
stateActivation = crocoddyl.ActivationModelWeightedQuad(np.matrix(stateWeights ** 2).T)
stateActivationTerm = crocoddyl.ActivationModelWeightedQuad(np.matrix(stateWeightsTerm ** 2).T)
xRegres = crocoddyl.ResidualModelState(state, defaultState)
xRegCost = crocoddyl.CostModelResidual(state, stateActivation, xRegres)
xRegTermCost = crocoddyl.CostModelResidual(state, stateActivationTerm, xRegres)

# input regulariztor
uResidual = crocoddyl.ResidualModelControl(state, actuation.nu)
uRegCost = crocoddyl.CostModelResidual(state, uResidual)

# Create cost model per each action model
runningCostModel = crocoddyl.CostModelSum(state, actuation.nu)
terminalCostModel = crocoddyl.CostModelSum(state, actuation.nu)

runningCostModel.addCost("gripperPose", goalTrackingCost, 1000)
runningCostModel.addCost("stateReg", xRegCost, 3e-1)
runningCostModel.addCost("ctrlReg", uRegCost, 1e-3)
runningCostModel.addCost("limitCost", limitCost, 1e3)
runningCostModel.addCost("obstacleAvoidance", obsAvoidCost, 2)

terminalCostModel.addCost("gripperPose", goalTrackingCost, 10000000)
terminalCostModel.addCost("stateReg", xRegTermCost, 1e-6)
terminalCostModel.addCost("limitCost", limitCost, 1e3)

runningCostModel2 = crocoddyl.CostModelSum(state, actuation.nu)
terminalCostModel2 = crocoddyl.CostModelSum(state, actuation.nu)

runningCostModel2.addCost("gripperPose", goalTrackingCost, 100)
runningCostModel2.addCost("stateReg", xRegCost, 1e-3)
runningCostModel2.addCost("ctrlReg", uRegCost, 1e-3)
runningCostModel2.addCost("limitCost", limitCost, 1e3)

terminalCostModel2.addCost("gripperPose", goalTrackingCost, 2000)
terminalCostModel2.addCost("stateReg", xRegTermCost, 1e-6)
terminalCostModel2.addCost("limitCost", limitCost, 1e3)

#%% Create Action function (dynamic model)
# Single integrator
runningModel=ActionModelRobot(state,actuation,runningCostModel)
terminalModel=ActionModelRobot(state,actuation,terminalCostModel)

runningModel.init_robot_sys(Single_integ(rmodel,dt))
terminalModel.init_robot_sys(Single_integ(rmodel,0))

runningModel2=ActionModelRobot(state,actuation,runningCostModel2)
terminalModel2=ActionModelRobot(state,actuation,terminalCostModel2)

runningModel2.init_robot_sys(Single_integ(rmodel,dt))
terminalModel2.init_robot_sys(Single_integ(rmodel,0))

# dynamic model
# dmodelrunning = crocoddyl.DifferentialActionModelFreeFwdDynamics(state, actuation, runningCostModel)
# dmodelTerminal = crocoddyl.DifferentialActionModelFreeFwdDynamics(state, actuation, terminalCostModel)

# runningModel = crocoddyl.IntegratedActionModelEuler(dmodelrunning, dt)
# terminalModel = crocoddyl.IntegratedActionModelEuler(dmodelTerminal, 0)

#%% Initialize the problem
problem0 = crocoddyl.ShootingProblem(x0, [runningModel2] * T, terminalModel2)
ddp0 = crocoddyl.SolverDDP(problem0)
ddp0.setCallbacks([crocoddyl.CallbackLogger(),crocoddyl.CallbackVerbose()])
ddp0.th_stop = 1e-3

problem = crocoddyl.ShootingProblem(x0, [runningModel] * T, terminalModel)
ddp = crocoddyl.SolverDDP(problem)
ddp.setCallbacks([crocoddyl.CallbackLogger(),crocoddyl.CallbackVerbose()])
ddp.th_stop = 1e-3

#%% Solve the problem
xs0 = [x0] * (T+1) # zero initialization
us0 =list(np.zeros((T, state.nv, 1))) # zero initialization

results = ddp0.solve(xs0,us0,1000)
if results == False:
    raise ValueError("ddp could not solve the shooting problem")
xs000 = ddp0.xs
xs00 = np.array(xs000).copy()

xs = ddp0.xs
us = ddp0.us

xs = [x0] * (T+1) # zero initialization
us = us0 =list(np.zeros((T, state.nv, 1)))

results = ddp.solve(xs,us,1000)
if results == False:
    raise ValueError("ddp could not solve the shooting problem")
    
#%% Start Pybullet
if p.getConnectionInfo()['isConnected'] == 1:
    p.disconnect() # disconnect from the server if it already have been started
p.connect(p.GUI)
p.configureDebugVisualizer(p.COV_ENABLE_GUI,0)
p.setAdditionalSearchPath(pybullet_data.getDataPath())
p.resetDebugVisualizerCamera(cameraDistance = camDistance, cameraTargetPosition = camTargetPos, cameraPitch = pitch, cameraYaw = yaw)
p.setTimeStep(dt)
p.loadURDF("plane.urdf")
p.setGravity(0,0,-9.81)
startPos = [0,0,0]
startOrientation = p.getQuaternionFromEuler([0,0,0])
robot_id = p.loadURDF(URDF, startPos, startOrientation, useFixedBase=1)#, flags = p.URDF_USE_SELF_COLLISION)
p.setRealTimeSimulation(0)

#%% Find the joint mapping between the pinocchio and pybullet (Unfertunately, they do not use the same numbers for the joints)
joints=[]
joints_name = []
pin_index = []
for i in range(p.getNumJoints(robot_id)):
 joint_info = np.array(p.getJointInfo(robot_id,i))
 if joint_info[2]==0 or joint_info[2]==1: # 0 for revolute and 1 for prismatic joints
     joints.append(i)
     joints_name.append(str(joint_info[1])[2:-1])
     pin_index.append(rmodel.getJointId(joints_name[-1]))
     
#%% Apply the solution on the robot
if movie_flag == 1:
    viewMatrix = p.computeViewMatrixFromYawPitchRoll(camTargetPos, camDistance, yaw, pitch, roll,
                                                         upAxisIndex)
    projectionMatrix = [
        1.0825318098068237, 0.0, 0.0, 0.0, 0.0, 1.732050895690918, 0.0, 0.0, 0.0, 0.0,
        -1.0002000331878662, -1.0, 0.0, 0.0, -0.020002000033855438, 0.0]
    out = cv.VideoWriter(movie_save, cv.VideoWriter_fourcc(*'DIVX'), int(1/dt), (pixelWidth,pixelHeight))
    
    
joint_limits = get_joint_limits(robot_id, joints) # to check if every thing is working fine
traj = np.array(ddp.xs)[:,:rmodel.nq]
traj0 = np.array(xs00)[:,:rmodel.nq]
goal_sphere_id=p.createVisualShape(p.GEOM_SPHERE,radius=0.05,rgbaColor=[0,1,0,1])
p.createMultiBody(baseVisualShapeIndex=goal_sphere_id,basePosition=[pos_goal[0],pos_goal[1],pos_goal[2]])

obs_sphere_id=p.createVisualShape(p.GEOM_SPHERE,radius=0.05,rgbaColor=[1,0,0,1])
p.createMultiBody(baseVisualShapeIndex=obs_sphere_id,basePosition=[obsPos[0],obsPos[1],obsPos[2]])

for i in range(obs_poses.shape[0]):
    rep = obs_poses[i,:]
    # obs_sphere_id=p.createVisualShape(p.GEOM_SPHERE,radius=0.03,rgbaColor=[1,0,0,1])
    #p.createMultiBody(baseVisualShapeIndex=p.createVisualShape(p.GEOM_SPHERE,radius=0.03,rgbaColor=[1,0,0,1]),basePosition=[rep[0],rep[1],rep[2]])
# projected objectives 
pos_goal2 = R @ pos_goal
print("original target", pos_goal, "projected target", pos_goal2)

pos_obs2 = R @ obsPos
print("original pos", obsPos, "projected obstacle ", pos_obs2)

for j in range(len(joints)):
        p.resetJointState(robot_id, joints[j], q0[rmodel.getJointId(joints_name[j])-1])
        p.setJointMotorControl2(robot_id, joints[j], controlMode = control_mode, force = max_force)


for t in range(len(xs)):
    q=traj[t,:]
    for j in range(len(joints)):
        p.setJointMotorControl2(robot_id, joints[j], controlMode = control_mode, targetPosition = q[rmodel.getJointId(joints_name[j]) - 1])
    p.stepSimulation()
    if movie_flag == 1:
        print(t)
        print("hellppppp")
        a = 2
        img_arr = p.getCameraImage(pixelWidth,
                               pixelHeight,pitch , 
                               viewMatrix=viewMatrix,
                               projectionMatrix=projectionMatrix,
                               shadow=1,
                               lightDirection=[1, 1, 1])
        img_bgr = img_arr[2][:,:,:3]
        img = cv.cvtColor(img_bgr,cv.COLOR_BGR2RGB)
        out.write(img)
    time.sleep(dt)

if movie_flag==1:
    out.release()   
    
    

#%% Plot user veiw

pos_proj, _ = pointProject(pos_goal.transpose(), yaw, pitch, camDistance, camTargetPos)
obs_proj, _ = pointProject(obsPos.transpose(), yaw, pitch, camDistance, camTargetPos)

plt.plot(pos_proj[0],-pos_proj[1],'go',markersize = 10)
plt.plot(obs_proj[0],-obs_proj[1],'ro',markersize = 10)

for i in range(obs_poses.shape[0]):
    rep = obs_poses[i,:]
    rep_proj, _ = pointProject(rep.transpose(), yaw, pitch, camDistance, camTargetPos)
    plt.plot(rep_proj[0],-rep_proj[1],'ro',markersize = 4)
ee_proj_x = []
ee_proj_y = []

ee0_proj_x = []
ee0_proj_y = []
for i in range(traj.shape[0]):
    q = traj[i,:]
    pin.forwardKinematics(rmodel, rdata, q, np.zeros([state.nv]))
    pin.updateFramePlacements(rmodel, rdata)
    ee_pos = rdata.oMf[eeIdx].translation
    ee_proj, _ = pointProject(ee_pos.transpose(), yaw, pitch, camDistance, camTargetPos)
    ee_proj_x.append(ee_proj[0]) 
    ee_proj_y.append(-ee_proj[1])
    
    qi = traj0[i,:]
    pin.forwardKinematics(rmodel, rdata, qi, np.zeros([state.nv]))
    pin.updateFramePlacements(rmodel, rdata)
    ee0_pos = rdata.oMf[eeIdx].translation
    ee0_proj, _ = pointProject(ee0_pos.transpose(), yaw, pitch, camDistance, camTargetPos)
    ee0_proj_x.append(ee0_proj[0]) 
    ee0_proj_y.append(-ee0_proj[1])
    
plt.plot(ee_proj_x,ee_proj_y)
plt.plot(ee0_proj_x,ee0_proj_y,'y',)
plt.xlim(-0.2,0.6)
plt.ylim(-0.2, 0.6)
plt.show()



pin.forwardKinematics(rmodel, rdata, q0, np.zeros([state.nv]))
pin.updateFramePlacements(rmodel, rdata)
fwd_x0T = rdata.oMf[eeIdx].translation
