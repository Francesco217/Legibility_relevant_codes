 function Manipulator_projection
% iLQR applied to a n-axis manipulator in 3D space with DH parameters

% addpath('./functions/');
% addpath('./projections/');

%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param.dt        = 1E-2; % Time step size
param.nbData    = 50;   % Number of datapoints
param.nbIter    = 100;  % Maximum number of iterations for iLQR
param.nbPoints  = 2;    % Number of viapoints
param.nbVarX    = 7;    % State space dimension (x1,x2,x3,...)
param.nbVarU    = param.nbVarX; % Control space dimension (dx1,dx2,dx3,...)
param.nbVarF2   = 6;    % Task space dimension (f1,f2,f3 for position, f4,f5,f6,f7 for unit quaternion)
param.nbVarF    = 3;
param.nbObj     = 2;    % number of possible targets
param.r         = 5;    % Control weighting term

Rtmp = q2R([cos(pi/3); [1;0;0]*sin(pi/3)]);
param.MuR = cat(3, Rtmp, Rtmp);

% Defining the point of view of the observer
param.pov    = [0.8;0;0];
param.ThetaX = 0;
param.ThetaY = 0;
param.ThetaZ = pi/2;

param.idx = [4,5,6, 8, 9];       % idx we care about (ignore basis and 3,7 bcs redundant)

param.Q1 = zeros(param.nbVarF*param.nbData, param.nbVarF*param.nbData, length(param.idx)-1);

% attractive weight is different depending on the joint
for k = 1:length(param.idx)-1
    Mat1 = eye(param.nbData);
    
    % terminal cost only for the end effector
    if k == length(param.idx)
        qf = param.nbData*1000;
        Mat1(param.nbData, param.nbData) = qf;
        param.Q1(:,:,k) = kron(Mat1, eye(param.nbVarF)*1000);
    else
    
        param.Q1(:,:,k) = kron(Mat1, eye(param.nbVarF)*10);
    end
end

% Weight matrix for the end effector
qf = param.nbData*100000;
Mat1 = eye(param.nbData);
Mat1(param.nbData, param.nbData) = qf;
param.Q2 = kron(Mat1, eye(param.nbVarF2)*100);

% Control matrix
param.R = speye((param.nbData-1)*param.nbVarU) * param.r;   


% DH parameters of Franka Emika robot 
param.dh.convention = 'm';                      % modified DH, a.k.a. Craig's formulation
param.dh.type = repmat('r', 1, param.nbVarX+1); % A rticulatory joints
param.dh.q = zeros(1, param.nbVarX+1);          % Angle about previous z
param.dh.q_offset = zeros(1, param.nbVarX+1);   % Offset on articulatory joints 
param.dh.alpha = [0 -pi/2, pi/2, pi/2, -pi/2, pi/2, pi/2, 0];   % Angle about common normal
param.dh.d = [0.333, 0, 0.316, 0, 0.384, 0, 0, 0.107];          % Offset along previous z to the common normal
param.dh.r = [0, 0, 0, 0.0825, -0.0825, 0, 0.088, 0];           % Length of the common normal 

% start position in joint space
q0 = [0; 0; 0; -pi/2; 0; pi/2; 0];


% start position in task space
[x03d,~] = fkin(q0, param);
x03d = x03d(1:3,:);             % only consider the x-y-z coordinates, ignore quaternions.
x0 = Point_Projection(x03d, param.ThetaX, param.ThetaY, param.ThetaZ, param.pov);

% objectives
param.Obj3d = [[0.2;0.5;0], [0.3; 0.2; 0.2]];   % objective leading to bad result
param.Obj3d = [[0.5;0.1;0.5], [0.5; -0.1;0.5]];  % objective leading to slighlty better result

% Projecting the objectives on viewspace
param.Obj = Point_Projection(param.Obj3d, param.ThetaX, param.ThetaY, param.ThetaZ, param.pov);

% computing the covariance matrix of the various goals
X_obj = param.Obj - mean(param.Obj, 2);

Cov = X_obj*X_obj'/(param.nbObj-1);
epsilon = 1E-4 * diag(ones(2));
Cov = Cov+epsilon;              % ensure not to have a singular matrix

% Compute the eigenvectos of the covariance matrix
[eig_vec, D] = eig(Cov);
eig_max = max(D, [], 'all');

% tuning the weights as a function of the eigenvalues
% start by sorting the eigenvalues and rearanging eigenvectors
[~,I] = sort(diag(D),'descend');
eig_vec = eig_vec(:,I);
D = diag([10,5])/eig_max;


% TUNING PARAMETERS ------------------------------------------------------

% Q3 matrix, keep away from unwanted trajectories
scaling_fact = 1;

param.w2 = linspace(0.05, 0.01, param.nbData);

q = zeros(1,param.nbData);
q(:) = scaling_fact;
param.Q3 = kron(diag(q), eig_vec*D*eig_vec');


% computing ideal goal position
param.qN = end_position_with_proj(param);
param.Obj1_orient = fkin(param.qN, param);

size_rec = 0.02;

% PLOT end position ------------------------------------------------------

h_init = figure;
hold on

% start by plotting the targets

draw_cube(param.Obj3d(:,1), size_rec, [0.4660 0.6740 0.1880]);

for i = 2:param.nbObj
    
   draw_cube(param.Obj3d(:,i), size_rec, [0.8500 0.3250 0.0980]); 
    
end

ftmp = fkin0(param.qN, param);
ftmp0 = fkin0(q0, param);

plot3(ftmp(1,:), ftmp(2,:), ftmp(3,:), '-','linewidth',4,'color',[.6 .6 .6]);
plot3(ftmp0(1,:), ftmp0(2,:), ftmp0(3,:), '-', 'linewidth', 4, 'color', [0 0 0]);
Rx = [1 0 0; 0 cos(param.ThetaX) -sin(param.ThetaX); 0 sin(param.ThetaX) cos(param.ThetaX)];
Ry = [cos(param.ThetaY) 0 sin(param.ThetaY); 0 1 0; -sin(param.ThetaY) 0 cos(param.ThetaY)];
Rz = [cos(param.ThetaZ) -sin(param.ThetaZ) 0; sin(param.ThetaZ) cos(param.ThetaZ) 0; 0 0 1];

ey = [0;0.1;0];
ey_new = Rz*Ry*Rx*ey;

quiver3(param.pov(1), param.pov(2), param.pov(3), ey_new(1), ey_new(2), ey_new(3),"LineWidth", 2, 'color', 'b');
draw_cube(param.pov, size_rec, [0.9290 0.6940 0.1250])

xlim([-0.5,1])
ylim([-0.5,1])
zlim([-0.5,1])
xlabel('x');
ylabel('y');
zlabel('z');
view(3)

hold off

% Plotting the 2D projection of the end position
h_init2 = figure;
hold on
rec1_x = [param.Obj(1,1)-size_rec, param.Obj(1,1)-size_rec,  param.Obj(1,1)+size_rec, param.Obj(1,1)+size_rec];
rec1_y = [param.Obj(2,1)-size_rec, param.Obj(2,1)+size_rec, param.Obj(2,1)+size_rec, param.Obj(2,1)-size_rec];

patch(rec1_x, rec1_y, [0.4660 0.6740 0.1880]);
for j = 2:param.nbObj
    rec2_x = [param.Obj(1,j)-size_rec, param.Obj(1,j)-size_rec, param.Obj(1,j)+size_rec, param.Obj(1,j)+size_rec];
    rec2_y = [param.Obj(2,j)-size_rec, param.Obj(2,j)+size_rec, param.Obj(2,j)+size_rec, param.Obj(2,j)-size_rec];
    patch(rec2_x, rec2_y, [0.8500 0.3250 0.0980]);
end

% plotting fastest path
plot([x0(1) param.Obj(1,2)], [x0(2) param.Obj(2,2)], 'color', [0.4940 0.1840 0.5560], "LineWidth", 0.8, "LineStyle", '--');
plot([x0(1) param.Obj(1,1)], [x0(2) param.Obj(2,1)], 'color', [0.4940 0.1840 0.5560], "LineWidth", 0.8, "LineStyle", '--');
f2d = Point_Projection(ftmp, param.ThetaX, param.ThetaY, param.ThetaZ, param.pov);
plot(f2d(1,:), f2d(2,:), '-','linewidth',4,'color',[.6 .6 .6]);
xlim([-0.5 1]);
ylim([-0.5, 1]);
hold off

% -------------------------------------------------------------------------

% defining the trajectories
param.traj3d = zeros(param.nbVarF, param.nbData, param.nbObj-1);
param.traj   = zeros(param.nbVarF-1, param.nbData, param.nbObj-1);

for k = 2:param.nbObj
   param.traj3d(1,:,k-1) = linspace(x03d(1), param.Obj3d(1,k), param.nbData);
   param.traj3d(2,:,k-1) = linspace(x03d(2), param.Obj3d(2,k), param.nbData);
   param.traj3d(3,:,k-1) = linspace(x03d(3), param.Obj3d(3,k), param.nbData);
   
   % Projecting the trajectory in the view space
   param.traj(:,:,k-1) = Point_Projection(param.traj3d(:,:,k-1), param.ThetaX, param.ThetaY, param.ThetaZ, param.pov);
end

% Time occurrence of viapoints
tl = linspace(1, param.nbData, param.nbPoints+1);
tl = round(tl(2:end));



%% Iterative LQR (iLQR)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u = zeros(param.nbVarU*(param.nbData-1), 1); %Initial control commands


% Transfer matrices (for linear system as single integrator)
Su0 = [zeros(param.nbVarX, param.nbVarX*(param.nbData-1)); kron(tril(ones(param.nbData-1)), eye(param.nbVarX)*param.dt)];
Sx0 = kron(ones(param.nbData,1), eye(param.nbVarX));

for n=1:param.nbIter	
    
	q = reshape(Su0 * u + Sx0 * q0, param.nbVarX, param.nbData); %System evolution
	
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
			break;1
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
Rx = [1 0 0; 0 cos(param.ThetaX) -sin(param.ThetaX); 0 sin(param.ThetaX) cos(param.ThetaX)];
Ry = [cos(param.ThetaY) 0 sin(param.ThetaY); 0 1 0; -sin(param.ThetaY) 0 cos(param.ThetaY)];
Rz = [cos(param.ThetaZ) -sin(param.ThetaZ) 0; sin(param.ThetaZ) cos(param.ThetaZ) 0; 0 0 1];

ey = [0;0.1;0];
ey_new = Rz*Ry*Rx*ey;

quiver3(param.pov(1), param.pov(2), param.pov(3), ey_new(1), ey_new(2), ey_new(3),"LineWidth", 2, 'color', 'b');
draw_cube(param.pov, size_rec, [0.9290 0.6940 0.1250])
xlim([-0.5,1]);
ylim([-0.5,1]);
zlim([-0.5,1]);

% Plot robot
pause(3)
for i = 1:param.nbData
    ftmp = fkin0(q(:,i), param);
    if i == 1
        hf = plot3(ftmp(1,:), ftmp(2,:), ftmp(3,:), '-','linewidth',4,'color',[.6 .6 .6]);
    else
        hf = plot3(ftmp(1,:), ftmp(2,:), ftmp(3,:), '-', 'linewidth', 4, 'color', [0.1, 0.1, 0.1]);
    end

    % Plot end-effector
    [ftmp, Rtmp] = fkin(q, param); 
    hf2 = plot3(ftmp(1,:), ftmp(2,:), ftmp(3,:), 'k-','linewidth',1);
    hf3 = plotCoordSys(ftmp(:,[1,tl]), Rtmp(:,:,[1,tl]) * .05);
    pause(0.1);
    if not(i == param.nbData || i==1)
        delete(hf)
        delete(hf2)
        delete(hf3)
    end
end

plotCoordSys(zeros(3,1), eye(3) * .1);
xlabel('f_1','fontsize',28); ylabel('f_2','fontsize',28); zlabel('f_3','fontsize',28);
xlim([-0.5,1]);
ylim([-0.5,1]);
zlim([-0.5,1]);

% Plot robot
view(3); axis vis3d; axis square; axis equal; 
[f_final,~] = fkin(q, param);
f_final = f_final(1:3,:);
h_2d = figure;
hold on 
x = Point_Projection(f_final, param.ThetaX, param.ThetaY, param.ThetaZ, param.pov);
rec1_x = [param.Obj(1,1)-size_rec, param.Obj(1,1)-size_rec,  param.Obj(1,1)+size_rec, param.Obj(1,1)+size_rec];
rec1_y = [param.Obj(2,1)-size_rec, param.Obj(2,1)+size_rec, param.Obj(2,1)+size_rec, param.Obj(2,1)-size_rec];

patch(rec1_x, rec1_y, [0.4660 0.6740 0.1880]);
for j = 2:param.nbObj
    rec2_x = [param.Obj(1,j)-size_rec, param.Obj(1,j)-size_rec, param.Obj(1,j)+size_rec, param.Obj(1,j)+size_rec];
    rec2_y = [param.Obj(2,j)-size_rec, param.Obj(2,j)+size_rec, param.Obj(2,j)+size_rec, param.Obj(2,j)-size_rec];
    patch(rec2_x, rec2_y, [0.8500 0.3250 0.0980]);
end
plot(x(1,:),x(2,:),"LineWidth", 1.5, "Marker", 'o')

% plotting fastest path
plot([x0(1) param.Obj(1,2)], [x0(2) param.Obj(2,2)], 'color', [0.4940 0.1840 0.5560], "LineWidth", 0.8, "LineStyle", '--');
plot([x0(1) param.Obj(1,1)], [x0(2) param.Obj(2,1)], 'color', [0.4940 0.1840 0.5560], "LineWidth", 0.8, "LineStyle", '--');



hold off


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
	[V,~] = eigs(K); %+eye(4)*1E-10
	q = [V(4,1); V(1:3,1)]; %for quaternions as [w,x,y,z]
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
	J = logmap(F2, F1) / e; %Error by considering manifold
end

%%%%%%%%%%%%%%%%%%%%%%
% Cost and gradient for a viapoints reaching task (in object coordinate system)
function [f, J] = f_reach(x, param)
	f = logmap(fkin(x, param), param.Obj1_orient); %Error by considering manifold
	J = []; 
	for t=1:size(x,2)
		J = blkdiag(J, Jkin_num(x(:,t), param));
	end
end

%%%%%%%%%%%%%%%%%%%%%%
% Logarithmic map for R^3 x S^3 manifold (with e in tangent space)
function e = logmap(f, f0)
	e = f(1:3,:) - f0(1:3,:); %Error on R^3
	for t=1:size(f,2)
		e(4:6,t) = logmap_S3(f(4:7,t), f0(4:7,t)); %Implementation with S3 (as in Martijn's thesis)
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

%%%%%%%%%%%%%%%%%%%%%%
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


    Nq  = param.nbVarX;
    Nx  = param.nbVarF;
    Nxp = param.nbVarF -1;
    Nf  = param.nbVarF2;
    
    cost = u'*param.R*u;
    
    for t = 1:param.nbData         % iterating over the points of the trajectory
        Q2 = param.Q2((t*Nf-Nf+1):t*Nf,t*Nf-Nf+1:t*Nf); 
        
        % computing the forward kinematics of the current configuration
        [F,~] = fkin0(q(:,t), param);
        
        % selecting only the articulations we care about
        F = F(:,param.idx);
        Fp = Point_Projection(F, param.ThetaX, param.ThetaY, param.ThetaZ, param.pov);
        
        [f, ~] = f_reach(q(:,t), param);
        
        cost = cost + f'*Q2*f;
        Q3 = param.Q3((t*Nxp-Nxp+1):t*Nxp,(t*Nxp-Nxp+1):t*Nxp);
        for k = 1:length(param.idx)-1
            Q1 = param.Q1((t*Nx-Nx+1):t*Nx,(t*Nx-Nx+1):t*Nx,k); 
        
            % attractive part of the cost for the current articulation
            
            X1 = F(:,k) - param.Obj3d(:,1);
            cost = cost + X1'*Q1*X1;
            
            for i = 1:param.nbObj-1        % iterating over points of the repulsive trajectory
               for j = 2:param.nbData-1      % iterating over the number of objectives

                        X2 = Fp(:,k) - param.traj(:,j,i);
                        cost = cost - param.w2(i)*log((X2'*Q3*X2));
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
    
    % Projection matrix
    Rx = [1 0 0; 0 cos(param.ThetaX) -sin(param.ThetaX); 0 sin(param.ThetaX) cos(param.ThetaX)];
    Ry = [cos(param.ThetaY) 0 sin(param.ThetaY); 0 1 0; -sin(param.ThetaY) 0 cos(param.ThetaY)];
    Rz = [cos(param.ThetaZ) -sin(param.ThetaZ) 0; sin(param.ThetaZ) cos(param.ThetaZ) 0; 0 0 1];
    P = inv(Rz*Ry*Rx);
    P = P([1,3],:);

    Nq  = param.nbVarX;
    Nx  = param.nbVarF;
    Nxp = param.nbVarF -1;
    Nf  = param.nbVarF2;
    
    % Jacobian dc/dq
    Jq = zeros(param.nbVarX, param.nbData);
    Jx = zeros(param.nbVarF, param.nbData);
    
    for t = 1:param.nbData
        
        Q2 = param.Q2((t*Nf-Nf+1):t*Nf,(t*Nf-Nf+1):t*Nf);
        
        % computing the forward kinematics for the end goal
        [f, Jf] = f_reach(q(:,t), param);
        
        Jq(:,t) = Jq(:,t) + Jf'*Q2*f;
        
        % computing the forward kinematics of all the articulations
        [F,~] = fkin0(q(:,t), param);
        F = F(:,param.idx);
        Fp = Point_Projection(F, param.ThetaX, param.ThetaY, param.ThetaZ, param.pov);
        
        J_all = Jkin_all(q(:,t), param);
        
        Q3 = param.Q3((t*Nxp-Nxp+1):t*Nxp,(t*Nxp-Nxp+1):t*Nxp);
        
        for k = 1:length(param.idx)-1
            Q1 = param.Q1((t*Nx-Nx+1):t*Nx,(t*Nx-Nx+1):t*Nx,k);
           
            X1 = F(:,k) - param.Obj3d(:,1);
            % Jacobian of the attractive part of the cost
            Jx(:,t) = Jx(:,t) + Q1*X1;
                      
           for j = 1:param.nbObj-1     % iterate over the number of different objectives
                for i = 2:param.nbData  % iterate over the trajectory 
              
                    X2 = Fp(:,k)-param.traj(:,i,j);
             
                    Jx(:,t) = Jx(:,t) - param.w2(i)*2*P'*Q3*X2/(X2'*Q3*X2);

                end
            end  
            
            Jq(:,t) = Jq(:,t) + J_all(:,:,k)'*Jx(:,t);
        end

    end
    Jq = reshape(Jq, param.nbVarX*param.nbData,1);

end

% HESSIAN COMPUTATION -----------------------------------------------------

function Hq = Hess(q, param)
%{ 
    Computes the hessian of the cost function with respect to q
    q: (NxM) trajectory in joint space
    Hq: (NMxNM) matrix containing the Hessian of each step
    
    
%}
   
    % Projection matrix
    Rx = [1 0 0; 0 cos(param.ThetaX) -sin(param.ThetaX); 0 sin(param.ThetaX) cos(param.ThetaX)];
    Ry = [cos(param.ThetaY) 0 sin(param.ThetaY); 0 1 0; -sin(param.ThetaY) 0 cos(param.ThetaY)];
    Rz = [cos(param.ThetaZ) -sin(param.ThetaZ) 0; sin(param.ThetaZ) cos(param.ThetaZ) 0; 0 0 1];
    P = inv(Rz*Ry*Rx);
    P = P([1,3],:);
    
    Nx  = param.nbVarF;
    Nxp = param.nbVarF -1;
    Nq  = param.nbVarX;
    Nf  = param.nbVarF2;
    
    % going from the Hessian of with respect to x to Hessian with respect
    % to q
    Hq = zeros(param.nbVarX*param.nbData);
    
    for t = 1:param.nbData
        J_all = Jkin_all(q(:,t), param);
        [~,Jf] = f_reach(q(:,t), param);
        Q3 = param.Q3((t*Nxp-Nxp+1):t*Nxp,(t*Nxp-Nxp+1):t*Nxp);
        
        % computing the forward kinematics of all the articulations
        [F,~] = fkin0(q(:,t), param);
        F = F(:,param.idx);
        Fp = Point_Projection(F, param.ThetaX, param.ThetaY, param.ThetaZ, param.pov);
        for k = 1:length(param.idx)-1 %iterating over the number of articulations 
            Hx_i = 2*param.Q1((t*Nx-Nx+1):t*Nx,(t*Nx-Nx+1):t*Nx,k);
            
            for j = 1:param.nbObj-1         % iterate over the number of different objectives
                for i = 2:param.nbData      % iterate over the trajectory    
            
                    
                    X2 = Fp(:,k)-param.traj(:,i,j);
                    Hx_i = Hx_i - param.w2(i)*P'*((2*Q3/((X2'*Q3*X2)) - 4*Q3*X2*(Q3*X2)'/((X2'*Q3*X2)^2)))*P; 
                end
            end
            
            Hq(t*Nq-(Nq-1):t*Nq,t*Nq-(Nq-1):t*Nq) = Hq(t*Nq-(Nq-1):t*Nq,t*Nq-(Nq-1):t*Nq) + J_all(:,:,k)'*Hx_i*J_all(:,:,k);        
        end
     
      
        Q2 = param.Q2((t*Nf-Nf+1):t*Nf,(t*Nf-Nf+1):t*Nf);
        Hq(t*Nq-(Nq-1):t*Nq,t*Nq-(Nq-1):t*Nq) = Hq(t*Nq-(Nq-1):t*Nq,t*Nq-(Nq-1):t*Nq) + Jf'*Q2*Jf;
    end


end


