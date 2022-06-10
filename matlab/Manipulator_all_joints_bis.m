function Manipulator_all_joints_bis

%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param.dt = 1E-2;                % Time step size
param.nbData = 50;              % Number of datapoints
param.nbIter = 200;             % Maximum number of iterations for iLQR
param.nbPoints = 2;             % Number of viapoints
param.nbVarX = 7;               % State space dimension (x1,x2,x3,...)
param.nbVarU = param.nbVarX;    % Control space dimension (dx1,dx2,dx3,...)
param.nbVarF = 3;
param.nbObj = 2;                % number of possible targets
param.r = 5;                    % Control weighting term
Rtmp = q2R([cos(pi/3); [1;0;0]*sin(pi/3)]);
param.MuR = cat(3, Rtmp, Rtmp);

param.idx = [6, 8, 9];       % idx we care about (ignore basis and 3,7 bcs redundant)



% defining the weight of the cost function

param.w2 = ones(param.nbData)*0.1; % trajectory avoidance weight
param.w2(1:30) = linspace(2,0.1,30);

param.Q1 = zeros(param.nbVarF*param.nbData, param.nbVarF*param.nbData, length(param.idx));

% attractive weight is different depending on the joint
for k = 1:length(param.idx)
    Mat1 = eye(param.nbData);
    
    % terminal cost only for the end effector
    if k == length(param.idx)
        qf = param.nbData*1000;
        Mat1(param.nbData, param.nbData) = qf;
        param.Q1(:,:,k) = kron(Mat1, eye(param.nbVarF)*100);
    else
    
        param.Q1(:,:,k) = kron(Mat1, eye(param.nbVarF)*0);
    end
end

param.R = speye((param.nbData-1)*param.nbVarU) * param.r;   %control matrix
 

% DH parameters of Franka Emika robot 
param.dh.convention = 'm';                      % modified DH, a.k.a. Craig's formulation
param.dh.type = repmat('r', 1, param.nbVarX+1); % Articulatory joints
param.dh.q = zeros(1, param.nbVarX+1);          % Angle about previous z
param.dh.q_offset = zeros(1, param.nbVarX+1);   % Offset on articulatory joints 
param.dh.alpha = [0 -pi/2, pi/2, pi/2, -pi/2, pi/2, pi/2, 0]; % Angle about common normal
param.dh.d = [0.333, 0, 0.316, 0, 0.384, 0, 0, 0.107]; % Offset along previous z to the common normal
param.dh.r = [0, 0, 0, 0.0825, -0.0825, 0, 0.088, 0]; % Length of the common normal 


% start position in joint space
q0 = [0; 0; 0; -pi/2; 0; pi/2; 0];

% start position in task space
[x0,~] = fkin(q0, param);
x0 = x0(1:3,:);             % only consider the x-y-z coordinates, ignore quaternions. 

% objectives
param.Obj3d = [[0.2;0.5;0], [0.3; 0.2; 0.2]]; 
param.Obj3d = [[0.5;0.1;0.5], [0.5; -0.1;0.5]];  

% defining the trajectories
param.traj3d = zeros(param.nbVarF, param.nbData, param.nbObj-1);
% computing ideal goal position
qN = end_position(param);

% visualize the ideal end position
figure('position',[10,10,600,600],'color',[1,1,1]);
hold on; 
view(3)

% start by plotting the targets
size_rec = 0.02;
draw_cube(param.Obj3d(:,1), size_rec, [0.4660 0.6740 0.1880]);

for i = 2:param.nbObj
    
   draw_cube(param.Obj3d(:,i), size_rec, [0.8500 0.3250 0.0980]); 
    
end

ftmp = fkin0(qN, param);

plot3(ftmp(1,:), ftmp(2,:), ftmp(3,:), '-','linewidth',4,'color',[.6 .6 .6]);
xlim([-0.5,1])
ylim([-0.5,1])
zlim([-0.5,1])
hold off;

for k = 2:param.nbObj
   param.traj3d(1,:,k-1) = linspace(x0(1), param.Obj3d(1,k), param.nbData);
   param.traj3d(2,:,k-1) = linspace(x0(2), param.Obj3d(2,k), param.nbData);
   param.traj3d(3,:,k-1) = linspace(x0(3), param.Obj3d(3,k), param.nbData);
end

% Time occurrence of viapoints
tl = linspace(1, param.nbData, param.nbPoints+1);
tl = round(tl(2:end));



%% Iterative LQR (iLQR)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u = zeros(param.nbVarU*(param.nbData-1), 1); % Initial control commands

% Transfer matrices (for linear system as single integrator)
Su0 = [zeros(param.nbVarX, param.nbVarX*(param.nbData-1)); kron(tril(ones(param.nbData-1)), eye(param.nbVarX)*param.dt)];
Sx0 = kron(ones(param.nbData,1), eye(param.nbVarX));

for n=1:param.nbIter	
    
	q = reshape(Su0 * u + Sx0 * q0, param.nbVarX, param.nbData); % System evolution
	
    Ju = 2*param.R*u;
    Jq = Gradient_q(q, param);
    
    Huu = 2*param.R;
    Hq = Hess(q, param);
    
	du = (Su0'*Hq*Su0 + Huu)^(-1)*(-Su0'*Jq-Ju);
	
	% Estimate step size with backtracking line search method
	alpha = 1;
	cost0 = cost_fct(q, u, param);
    
    
    
	while 1
		utmp = u + du * alpha;
		qtmp = reshape(Su0 * utmp + Sx0 * q0, param.nbVarX, param.nbData);
        
		cost = cost_fct(qtmp, utmp, param);
		if cost < cost0 || alpha < 1E-3
			break;
		end
		alpha = alpha * 0.5;
	end
	u = u + du * alpha;
    clc;
    display(n);
	
	if norm(du * alpha) < 1E-2
		break; % Stop iLQR when solution is reached
	end
end
disp(['iLQR converged in ' num2str(n) ' iterations.']);

[x3d,~] =fkin(q(:,param.nbData), param); 
if abs(x3d(1) - param.Obj3d(1,1)) < 1E-5
    if abs(x3d(2) -param.Obj3d(2,1)) < 1E-5
        if abs(x3d(3)-param.Obj3d(3,1)) < 1E-5
        
            display("reached goal point");
        end
    end
end

%% Plot state space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = figure('position',[10,10,600,600],'color',[1,1,1]); hold on; rotate3d on; %axis off;

% start by plotting the targets
size_rec = 0.02;
draw_cube(param.Obj3d(:,1), size_rec, [0.4660 0.6740 0.1880]);

for i = 2:param.nbObj
    
   draw_cube(param.Obj3d(:,i), size_rec, [0.8500 0.3250 0.0980]); 
    
end
xlim([-0.5,1]);
ylim([-0.5,1]);
zlim([-0.5,1]);

% Plot robot
pause(3)
for i = 1:param.nbData
    ftmp = fkin0(q(:,i), param);
    if i == 1
        hf = plot3(ftmp(1,:), ftmp(2,:), ftmp(3,:), '-','linewidth',4,'color',[.6 .6 .6]);
    elseif i == 20
        hf = plot3(ftmp(1,:), ftmp(2,:), ftmp(3,:), '-', 'linewidth', 4, 'color', [0.3, 0.3, 0.3]);
    else
        hf = plot3(ftmp(1,:), ftmp(2,:), ftmp(3,:), '-', 'linewidth', 4, 'color', [0.1, 0.1, 0.1]);
    end
       
    % Plot end-effector
    [ftmp, Rtmp] = fkin(q, param); 
    hf2 = plot3(ftmp(1,:), ftmp(2,:), ftmp(3,:), 'k-','linewidth',1);
    hf3 = plotCoordSys(ftmp(:,[1,tl]), Rtmp(:,:,[1,tl]) * .05);
    pause(0.1);
    if not(i == param.nbData || i==1 || i == 20)
        delete(hf)
        delete(hf2)
        delete(hf3)
    end
end

% plot3(0, 0, 0, 'k.','markersize',50);
plotCoordSys(zeros(3,1), eye(3) * .1);
xlabel('f_1','fontsize',28); ylabel('f_2','fontsize',28); zlabel('f_3','fontsize',28);
xlim([-0.5,1]);
ylim([-0.5,1]);
zlim([-0.5,1]);
% Plot robot
view(3); axis vis3d; axis square; axis equal; 


end 

%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot coordinate system
function hf = plotCoordSys(x, R)
	hf = [];
	for t=1:size(x,2)
		hf = [hf, plot3([x(1,t), x(1,t)+R(1,1,t)], [x(2,t), x(2,t)+R(2,1,t)], [x(3,t), x(3,t)+R(3,1,t)], 'r-','linewidth',4)];
		hf = [hf, plot3([x(1,t), x(1,t)+R(1,2,t)], [x(2,t), x(2,t)+R(2,2,t)], [x(3,t), x(3,t)+R(3,2,t)], 'g-','linewidth',4)];
		hf = [hf, plot3([x(1,t), x(1,t)+R(1,3,t)], [x(2,t), x(2,t)+R(2,3,t)], [x(3,t), x(3,t)+R(3,3,t)], 'b-','linewidth',4)];
		hf = [hf, plot3(x(1,t), x(2,t), x(3,t), 'k.','markersize',20)];
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% Unit quaternion to rotation matrix conversion (for quaternions as [w,x,y,z])
function R = q2R(q)
	R = [1 - 2 * q(3)^2 - 2 * q(4)^2, 2 * q(2) * q(3) - 2 * q(4) * q(1), 2 * q(2) * q(4) + 2 * q(3) * q(1);
		 2 * q(2) * q(3) + 2 * q(4) * q(1), 1 - 2 * q(2)^2 - 2 * q(4)^2,  2 * q(3) * q(4) - 2 * q(2) * q(1);
		 2 * q(2) * q(4) - 2 * q(3) * q(1), 2 * q(3) * q(4) + 2 * q(2) * q(1), 1 - 2 * q(2)^2 - 2 * q(3)^2];
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% Rotation matrix to unit quaternion conversion 
function q = R2q(R)
	R = R';
	K = [R(1,1)-R(2,2)-R(3,3), R(2,1)+R(1,2), R(3,1)+R(1,3), R(2,3)-R(3,2);
		 R(2,1)+R(1,2), R(2,2)-R(1,1)-R(3,3), R(3,2)+R(2,3), R(3,1)-R(1,3);
		 R(3,1)+R(1,3), R(3,2)+R(2,3), R(3,3)-R(1,1)-R(2,2), R(1,2)-R(2,1);
		 R(2,3)-R(3,2), R(3,1)-R(1,3), R(1,2)-R(2,1), R(1,1)+R(2,2)+R(3,3)] / 3; 
	[V,~] = eigs(K); 
	q = [V(4,1); V(1:3,1)]; % for quaternions as [w,x,y,z]
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward kinematics for all robot articulations, for articulatory joints and 
% modified DH convention (a.k.a. Craig's formulation)
function [f, R] = fkin0(x, param)
	N = length(param.dh.q);
	x = [x; 0];
	Tf = eye(4);
	R = zeros(3,3,N);
	f = zeros(3,1);
	for n=1:N
		ct = cos(x(n) + param.dh.q_offset(n));
		st = sin(x(n) + param.dh.q_offset(n));
		ca = cos(param.dh.alpha(n));
		sa = sin(param.dh.alpha(n));
		Tf = Tf * [ct,    -st,     0,   param.dh.r(n)   ; ...
				   st*ca,  ct*ca, -sa, -param.dh.d(n)*sa; ...
				   st*sa,  ct*sa,  ca,  param.dh.d(n)*ca; ...
				   0,      0,      0,   1               ];
		R(:,:,n) = Tf(1:3,1:3); 
		f =  [f, Tf(1:3,end)];
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward kinematics for end-effector (from DH parameters)
function [f, R] = fkin(x, param)
	f = zeros(7, size(x,2));
	R = zeros(3, 3, size(x,2));
	for t=1:size(x,2)
		[ftmp, Rtmp] = fkin0(x(:,t), param);
		R(:,:,t) = Rtmp(:,:,end); 
		f(:,t) = [ftmp(:,end); R2q(Rtmp(:,:,end))];
	end
end

%%%%%%%%%%%%%%%%%%%%%%
% Jacobian of forward kinematics function with numerical computation
function J = Jkin_num(x, param)
	e = 1E-6;
	X = repmat(x, [1, param.nbVarX]);
	F1 = fkin(X, param);
	F2 = fkin(X + eye(param.nbVarX) * e, param);
	J = logmap(F2, F1) / e; % Error by considering manifold
end

%%%%%%%%%%%%%%%%%%%%%%
% Cost and gradient for a viapoints reaching task (in object coordinate system)
function [f, J] = f_reach(x, param)
	f = logmap(fkin(x, param), param.Mu); % Error by considering manifold
	J = []; 
	for t=1:size(x,2)
		J = blkdiag(J, Jkin_num(x(:,t), param));
	end
end

%%%%%%%%%%%%%%%%%%%%%%
% Logarithmic map for R^3 x S^3 manifold (with e in tangent space)
function e = logmap(f, f0)
	e = f(1:3,:) - f0(1:3,:); 
	for t=1:size(f,2)
		e(4:6,t) = logmap_S3(f(4:7,t), f0(4:7,t)); % Implementation with S3 (as in Martijn's thesis)
	end
end

%%%%%%%%%%%%%%%%%%%%%%
% Logarithmic map for S^3 manifold (with e in tangent space)
function u = logmap_S3(x, x0)
	R = QuatMatrix(x0);
	x = R' * x;
 	sc = acoslog(x(1,:)) ./ sqrt(1 - x(1,:).^2);
 	
 	sc(isnan(sc)) = 1;
	u = x(2:end,:) .* sc;
end

%%%%%%%%%%%%%%%%%%%%%%
% Matrix form of quaternion
function Q = QuatMatrix(q)
	Q = [q(1) -q(2) -q(3) -q(4);
	     q(2)  q(1) -q(4)  q(3);
	     q(3)  q(4)  q(1) -q(2);
	     q(4) -q(3)  q(2)  q(1)];
end


% Arcosine redefinition to make sure the distance between antipodal quaternions is zero (see Section 2.50 from Dubbelman's Thesis)
function y = acoslog(x)
	y = acos(x);
	id = (x>=-1.0 && x<0);
	y(id) = y(id) - pi;
end




%% Functions for Sd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Logarithmic map for S^d manifold (with e in tangent space)
function e = logmap_Sd(x, x0)
	for t=1:size(x,2)
		etmp = logmap2(x(:,t), x0(:,t));
		
		R = coordTangentSpace(x0(:,t));
		e(:,t) = R(1:end-1,:) * etmp;	
	end
end

% Logarithmic map for S^d manifold (with e in ambient space)
function e = logmap2(x, x0)
	for t=1:size(x,2)

		th = acoslog(x0(:,t)' * x(:,t));
		u = x(:,t) - x0(:,t)' * x(:,t) * x0(:,t);
		u = th .* u ./ norm(u);

		e(:,t) = u;
	end
end

% Parallel transport (see https://ronnybergmann.net/mvirt/manifolds/Sn.html, https://towardsdatascience.com/geodesic-regression-d0334de2d9d8)
function v2 = transp2(x, y, v)
	d = acoslog(x' * y);
	
	v2 = v - (logmap2(y,x) + logmap2(x,y)) .* (logmap2(y,x)' * v) ./ d.^2;
end

% Coordinate system of tangent space 
function R = coordTangentSpace(x)
	nbVar = size(x,1);
	e0 = [0; 0; 0; 1]; % Origin of the manifold
	R = zeros(nbVar);
	for n=1:nbVar-1
		v = zeros(nbVar,1);
		v(n) = 1;
		R(n,:) = transp2(e0, x, v);
	end
	R(nbVar,:) = x;
end

% Jacobian of forward kinematics for all articulation

function J = Jkin_all(q,param)

%{
    q:  (nbVarX x 1) column vector of the joint space coordinates
    param:  structure containing the parameters of the problem

    J: (nbVarF x nbVarX x nbVarX-1) matrix of matrix where the first
    two dimention are the Jacobian for one of the articulations. We fon't
    compute the forward of the first articulation as it's position is fixed

%}




e = 1E-6;
Q = repmat(q, [1, param.nbVarX]);
Qe = Q + eye(param.nbVarX) * e;

J = zeros(param.nbVarF, param.nbVarX, length(param.idx));


for i = 1:param.nbVarX
    
    [F1, ~] = fkin0(Q(:,i), param);
    [F2, ~] = fkin0(Qe(:,i), param);  
    
    for k = 1:length(param.idx)
    
        J(:,i,k) = (F2(:,param.idx(k))- F1(:,param.idx(k))) / e;
        
    end        
end

end

% COST FUNCTION -----------------------------------------------------------

function cost = cost_fct(q, u, param)
%{
    q:      (nbVarX x nbData) trajectory in joint space
    u:      (nbVarU x nbData-1) commands of the trajectory
    param:  structure containing the parameters of the problem
    cost:   (scalar) cost of trajectory
%}


   
    Nx = param.nbVarF;
    
    
    cost = u'*param.R*u;
    
    for t = 1:param.nbData         % iterating over the points of the trajectory
        % computing the forward kinematics of the current 
        [F,~] = fkin0(q(:,t), param);
        % selecting only the articulations we care about
        F = F(:,param.idx);
        
        for k = 1:length(param.idx)
            Q1 = param.Q1((t*Nx-Nx+1):t*Nx,(t*Nx-Nx+1):t*Nx,k); 
        
            % attractive part of the cost for the current articulation
            
            X1 = F(:,k) - param.Obj3d(:,1);
            cost = cost + X1'*Q1*X1;
            
            for i = 1:param.nbObj-1        % iterating over points of the repulsive trajectory
               for j = 2:param.nbData-1      % iterating over the number of objectives

                        X2 = F(:,k) - param.traj3d(:,j,i);
                        cost = cost - param.w2(i)*sum(abs(X2));
                end

            end
        end
    end

end

% GRADIENT COMPUTATION OF COST FUNCTION
function Jq = Gradient_q(q, param)
    %{ 
    Computes the gradient of the cost function with respect to q
    q: (NxM) trajectory in joint space
    Jx: (NMx1) vector containing the gradient of each step

    %}
    
     
    Nx = param.nbVarF;
    % Jacobian dc/dq
    Jq = zeros(param.nbVarX, param.nbData);
    Jx = zeros(param.nbVarF, param.nbData);
    
    
    for t = 1:param.nbData
        % computing the forward kinematics of all the articulations
        [F,~] = fkin0(q(:,t), param);
        F = F(:,param.idx);
        
        J_all = Jkin_all(q(:,t), param);
        
        for k = 1:length(param.idx)
            Q1 = param.Q1((t*Nx-Nx+1):t*Nx,(t*Nx-Nx+1):t*Nx,k);
           
            X1 = F(:,k) - param.Obj3d(:,1);
            % Jacobian of the attractive part of the cost
            Jx(:,t) = Jx(:,t) + Q1*X1;
                      
           for j = 1:param.nbObj-1     % iterate over the number of different objectives
                for i = 2:param.nbData  % iterate over the trajectory 

                    X2 = F(:,k)-param.traj3d(:,i,j);
                    I_neg = X2 <= 0;
                    I_pos = X2 > 0;
                    Jx(I_neg,t) = Jx(I_neg,t) + param.w2(i);
                    Jx(I_pos,t) = Jx(I_pos,t) - param.w2(i);

                end
            end  
            
            Jq(:,t) = Jq(:,t) + J_all(:,:,k)'*Jx(:,t);
        end
    end
    Jq = reshape(Jq, param.nbVarX*param.nbData,1);

end

% HESSIAN COMPUTATION

function Hq = Hess(q, param)
%{ 
    Computes the hessian of the cost function with respect to q
    q: (NxM) trajectory in joint space
    Hq: (NMxNM) matrix containing the Hessian of each step
    
    
%}
   
    
    Nx = param.nbVarF;
    Nq = param.nbVarX;
    
    %going from the Hessian of with respect to x to Hessian with respect to
    %q
    Hq = zeros(param.nbVarX*param.nbData);
    
    for i = 1:param.nbData
        J_all = Jkin_all(q(:,i), param);
        
        for k = 1:length(param.idx)
                        
            Hx_i = 2*param.Q1((i*Nx-Nx+1):i*Nx,(i*Nx-Nx+1):i*Nx,k);
            % computing the Jacobian of the kinematics for the current point
            Hq(i*Nq-(Nq-1):i*Nq,i*Nq-(Nq-1):i*Nq) = Hq(i*Nq-(Nq-1):i*Nq,i*Nq-(Nq-1):i*Nq) + J_all(:,:,k)'*Hx_i*J_all(:,:,k);
            
        end
    end


end


