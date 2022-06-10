    %%    iLQR applied to a planar manipulator for a viapoints task (batch formulation)
%%
%%    Copyright (c) 2021 Idiap Research Institute, http://www.idiap.ch/
%%    Written by Sylvain Calinon <https://calinon.ch>
%%
%%    This file is part of RCFS.
%%
%%    RCFS is free software: you can redistribute it and/or modify
%%    it under the terms of the GNU General Public License version 3 as
%%    published by the Free Software Foundation.
%%
%%    RCFS is distributed in the hope that it will be useful,
%%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%%    GNU General Public License for more details.
%%
%%    You should have received a copy of the GNU General Public License
%%    along with RCFS. If not, see <http://www.gnu.org/licenses/>.

function trajectory_avoidance_decrease_weight

% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param.dt = 1E-2;        %Time step length
param.nbData = 50;      %Number of datapoints
param.nbIter = 3000;    %Number of iterations for iLQR
param.nbPoints = 1;     %Number of viapoints
param.nbVarX = 3;       %State space dimension (x1,x2,x3)
param.nbVarU = 3;       %Control space dimension (dx1,dx2,dx3)
param.nbVarF = 2;       %Objective function dimension (f1,f2,f3, with f3 as orientation)
param.l = [2; 2; 1];    %Robot links lengths
param.sz = [.2, .3];    %Size of objects
param.r = 100;          %Control weight term
param.nbObj = 2;

% defining the two objectives
param.Obj3d = [[20.; -1.;0], [20.; 1.;0]];%, [20;0;-1]];
%param.Obj3d = [[22.; -4.], [5.; 5.]];
%param.Obj3d = [[15.; 6.], [-2.; 3.]];
%param.Obj3d = [[24; 7;5],[23; 4;10]];

%Initial robot pose
x0_3d = [0;0;0];


% generate the trajectory to be avoided for all objects exept goal (index 1)
param.traj3d = zeros(param.nbVarX, param.nbData, param.nbObj-1);

for i = 2:param.nbObj 
    param.traj3d(1,:,i-1) = linspace(0,param.Obj3d(1,i),param.nbData);
    param.traj3d(2,:,i-1) = linspace(0,param.Obj3d(2,i),param.nbData);
    param.traj3d(3,:,i-1) = linspace(0,param.Obj3d(3,i),param.nbData);
end


% -------------------------------------------------------------------------
% PLOTTING THE INITIAL SITUATION
figure
hold on
scatter3(x0_3d(1),x0_3d(2),x0_3d(3));
scatter3(param.Obj3d(1,:),param.Obj3d(2,:), param.Obj3d(3,:));

for j = 2:param.nbObj
    scatter3(param.traj3d(1,:,j-1), param.traj3d(2,:,j-1), param.traj3d(3,:,j-1));
    plot3([x0_3d(1) param.Obj3d(1,j)],[x0_3d(2) param.Obj3d(2,j)],[x0_3d(3) param.Obj3d(3,j)]);
end

size_rec = 0.5;

draw_cube(param.Obj3d(:,1), size_rec, 'green');

for i = 2:param.nbObj
    
   draw_cube(param.Obj3d(:,i), size_rec, 'red'); 
    
end

xlabel('x');
ylabel('y');
zlabel('z');
xlim([-2,30]);
ylim([-2,30]);
zlim([-2,30]);
title('initial situation');
grid on;
view(3);
hold off

% ------------------------------------------------------------------------

% defining the point of view of the observer
param.pov = [0; 0; 30];
param.ThetaX = -pi/2;
param.ThetaY = 0;
param.ThetaZ = 0;

% projecting objectives and trajectories on viewplane
param.Obj = Point_Projection(param.Obj3d,param.ThetaX, param.ThetaY, param.ThetaZ, param.pov);
param.traj = zeros(param.nbVarX-1, param.nbData, param.nbObj-1);

for i = 1:param.nbObj-1
    param.traj(:,:,i) = Point_Projection(param.traj3d(:,:,i), param.ThetaX, param.ThetaY, param.ThetaZ, param.pov);
end

% computing the covariance matrix of the various goals
X_obj = param.Obj - mean(param.Obj, 2);

Cov = X_obj*X_obj'/(param.nbObj-1);
epsilon = 1E-4 * diag(ones(2));
Cov = Cov+epsilon; % ensure not to have a singular matrixx

% Compute the eigenvectos of the covariance matrix
[eig_vec, D] = eig(Cov);
eig_max = max(D, [], 'all');

% tuning the weights as a function of the eigenvalues
% start by sorting the eigenvalues and rearanging eigenvectors
[~,I] = sort(diag(D),'descend');
eig_vec = eig_vec(:,I);
D = diag([13,1])/eig_max;



% TUNING PARAMETERS ------------------------------------------------------

% Q2 matrix, keep away from unwanted trajectories
scaling_fact = 30;
param.w2 = linspace(55, 1, param.nbData);

q = zeros(1,param.nbData);
q(:) = param.w2;
param.Q2 = kron(diag(q), eig_vec*D*eig_vec');


% Q1 brings end effector to right goal
param.Qf = param.nbData*param.w2(1)*scaling_fact;
Mat1 = eye(param.nbData);
Mat1(param.nbData, param.nbData) = param.Qf; %terminal weight
param.Q1 = kron(Mat1, eye(param.nbVarX)*10000);

% Qp brings end effector to right goal in projection space
Mat2 = eye(param.nbData);
Mat2(param.nbData, param.nbData) = param.Qf;
param.Qp = kron(Mat2, eig_vec*eye(param.nbVarX-1)*eig_vec');


% Defining the control matrix
R = speye((param.nbData-1)*param.nbVarU)*param.r;

%--------------------------------------------------------------------------

% Iterative LQR (iLQR)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u = zeros(param.nbVarU*(param.nbData-1), 1); %Initial commands

% projecting it on viewspace
x0 = Point_Projection(x0_3d, param.ThetaX, param.ThetaY, param.ThetaZ, param.pov);

% Transfer matrices (for linear system as single integrator)
Su0 = [zeros(param.nbVarX, (param.nbVarX)*(param.nbData-1)); tril(kron(ones(param.nbData-1), eye(param.nbVarX)*param.dt))];
Sx0 = kron(ones(param.nbData,1), eye(param.nbVarX));

for n=1:param.nbIter	
	x3d = reshape(Su0 * u + Sx0 * x0_3d, param.nbVarX, param.nbData); %System evolution  
    
    % computing the Jacobians and Hessian of the cost function around the
    % current trajectory
    Jx = Gradient_x(x3d, param);
    Ju = 2*R*u;
    
    Hxx = Hess_x(x3d, param);
    Huu = 2*R;
    
    % calculating the update step
    du = (Su0'*Hxx*Su0 + Huu)^(-1)*(-Su0'*Jx-Ju);
    
	% Estimate step size with backtracking line search method
	
	alpha = 1;
	cost0 = cost_fct(x3d, u, param, R);
    
	while 1
		utmp = u + du * alpha;
		xtmp = reshape(Su0 * utmp + Sx0 * x0_3d, param.nbVarX, param.nbData); % System evolution
        
        cost = cost_fct(xtmp, utmp, param, R);
        
        if cost < cost0 || alpha < 1E-3
			break;
        end
        
		alpha = alpha * 0.5;
    end
    
	u = u + du * alpha;
    
    clc;
    display(n);
	
	if norm(du * alpha) < 1E-3
		break; %Stop iLQR when solution is reached
    end
    
end
disp(['iLQR converged in ' num2str(n) ' iterations.']);


%% Plot state space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if abs(x3d(1,param.nbData) - param.Obj(1,1)) < 1E-5
    if abs(x3d(2, param.nbData) -param.Obj(2,1)) < 1E-5
        if abs(x3d(2, param.nbData)-param.Obj(2,1)) < 1E-5
        
            display("reached goal point");
        end
    end
end

x = Point_Projection(x3d, param.ThetaX, param.ThetaY, param.ThetaZ, param.pov);
figure('position',[10,10,550,550], 'color', [1 1 1]);
axis off;
hold on;

Xs = linspace(-50,50,200);
[Grid_X, Grid_Y] = meshgrid(Xs);
Grid_X = reshape(Grid_X, 1, []);
Grid_Y = reshape(Grid_Y, 1, []);
Grid_coord = [Grid_X; Grid_Y];
costP = cost_plot(Grid_coord, param);
costP = reshape(costP, 200, 200);

size_rec = 0.5;
rec1_x = [param.Obj(1,1)-size_rec, param.Obj(1,1)-size_rec,  param.Obj(1,1)+size_rec, param.Obj(1,1)+size_rec];
rec1_y = [param.Obj(2,1)-size_rec, param.Obj(2,1)+size_rec, param.Obj(2,1)+size_rec, param.Obj(2,1)-size_rec];

patch(rec1_x, rec1_y, [0.4660 0.6740 0.1880]);
for j = 2:param.nbObj
    rec2_x = [param.Obj(1,j)-size_rec, param.Obj(1,j)-size_rec, param.Obj(1,j)+size_rec, param.Obj(1,j)+size_rec];
    rec2_y = [param.Obj(2,j)-size_rec, param.Obj(2,j)+size_rec, param.Obj(2,j)+size_rec, param.Obj(2,j)-size_rec];
    patch(rec2_x, rec2_y, [0.8500 0.3250 0.0980]);
end
plot(x(1,:),x(2,:), "LineWidth", 1.5, "Marker", 'o')

% plotting fastest path
plot([x0(1) param.Obj(1,2)], [x0(2) param.Obj(2,2)], 'color', [0.4940 0.1840 0.5560], "LineWidth", 0.8, "LineStyle", '--');
plot([x0(1) param.Obj(1,1)], [x0(2) param.Obj(2,1)], 'color', [0.4940 0.1840 0.5560], "LineWidth", 0.8, "LineStyle", '--');
xlim([-1,22])
ylim([-12.5,12.5])
legend('Location', 'best')
legend("desired target", "undesired target", "legible trajectory", "optimal trajectory");

xlabel('x');
ylabel('y');

grid on
hold off;

figure
hold on
grid on

plot3(x3d(1,:),x3d(2,:),x3d(3,:));
scatter3(x0_3d(1), x0_3d(2), x0_3d(3));
scatter3(param.Obj3d(1,:), param.Obj3d(2,:), param.Obj3d(3,:));

% plotting the objectives
draw_cube(param.Obj3d(:,1), size_rec, 'green');

for i = 2:param.nbObj
    
   draw_cube(param.Obj3d(:,i), size_rec, 'red'); 
    
end



xlim([-2,30]);
ylim([-2,30]);
zlim([-2,30]);
xlabel('x');
ylabel('y');
zlabel('z');

view(3);
hold off;
figure
Xs = linspace(-50,50,200);
[Grid_X, Grid_Y] = meshgrid(Xs);
mesh(Grid_X, Grid_Y,costP)

end



function Jx = Gradient_x(X, param)
    %{ 
    Computes the gradient of the cost function with respect to the x
    X: (NxM) datapoints
    Jx: (NMx1) vector containing the gradient of each step

    %}
    
    % building projection matrix
    Rx = [1 0 0; 0 cos(param.ThetaX) -sin(param.ThetaX); 0 sin(param.ThetaX) cos(param.ThetaX)];
    Ry = [cos(param.ThetaY) 0 sin(param.ThetaY); 0 1 0; -sin(param.ThetaY) 0 cos(param.ThetaY)];
    Rz = [cos(param.ThetaZ) -sin(param.ThetaZ) 0; sin(param.ThetaZ) cos(param.ThetaZ) 0; 0 0 1];
    P = inv(Rz*Ry*Rx);
    P = P([1,3],:);
    
    X1 = reshape(X-param.Obj3d(:,1), param.nbVarX*param.nbData,1);
    Q2 = param.Q2;
    Jx = 2*param.Q1*X1;
    
    % Projecting X on 2D viewplane to calcultate 2D Jacobian in viewplane
    Xp = Point_Projection(X, param.ThetaX, param.ThetaY, param.ThetaZ, param.pov);
    X1p = reshape(Xp-param.Obj(:,1), (param.nbVarX-1)*param.nbData, 1);
    
    Jxp = 2*param.Qp*X1p;
    
    for j = 1:param.nbObj-1     % iterate over the number of different objectives
        for i = 2:param.nbData  % iterate over the trajectory 
                
            X2 = reshape(Xp-param.traj(:,i,j), (param.nbVarX-1)*param.nbData,1);
            Jxp = Jxp - 2*param.Q2*X2;
            
        end
    end   
    
    % going from 2d Jacobian to the 3D one
    Jxp = reshape(Jxp, param.nbVarX-1, param.nbData);
    Jxp = reshape(P'*Jxp, param.nbVarX*param.nbData,1);
    Jx = Jx + Jxp;
 
end

function Jx = Gradient_plot(X,param)
    
    % -- NOT UP TO DATE: must update gradient calculation accountinf for
    % projection 

    Q1 = param.Q1(1:2,1:2);
    Q2 = param.Q2(1:2,1:2);
    [N,M] = size(X);
    Jx = zeros(N,M);
    
    for k = 1:M
       X1 = X(:,k) - param.Obj(:,1);
       Jx(:,k) =  Q1*X1;
       for i = 1:50
           
           for j = 1:param.nbObj-1
               X2 = X(:,k) - param.traj(:,i,j);
               Jx(:,k) = Jx(:,k)-param.w2(i)*2*Q2*X2/(X2'*Q2*X2);
           end
       end
        
        
    end

end

function Hxx = Hess_x(X, param)
    %{
    Returns the Hessian of the cost function around each one of the
    datapoints
    X (NxM) input data
    Hxx (N*NbrDataPoints) square matrix
    %}
    
    % building projection matrix
    Rx = [1 0 0; 0 cos(param.ThetaX) -sin(param.ThetaX); 0 sin(param.ThetaX) cos(param.ThetaX)];
    Ry = [cos(param.ThetaY) 0 sin(param.ThetaY); 0 1 0; -sin(param.ThetaY) 0 cos(param.ThetaY)];
    Rz = [cos(param.ThetaZ) -sin(param.ThetaZ) 0; sin(param.ThetaZ) cos(param.ThetaZ) 0; 0 0 1];
    P = inv(Rz*Ry*Rx);
    P = P([1,3],:);
      
    Q1 = param.Q1;

    Hxx = 2*Q1;
    
    % projecting trajectory on the viewplane
    Xp = Point_Projection(X, param.ThetaX, param.ThetaY, param.ThetaZ, param.pov);
    Hxxp = []; % Hessian for all points 
    
    for k = 1:param.nbData
        H2d = 2*param.Qp(2*k-1:2*k, 2*k-1:2*k);
        
        for j = 1:param.nbObj-1         % iterate over the number of different objectives
            for i = 2:param.nbData      % iterate over the trajectory    
            
                % computing Hessian in projection space
                Q2 = param.Q2(2*k-1:2*k,2*k-1:2*k);
                H2d = H2d - 2*Q2; 
                
            end
        end
        
        H3d = P'*H2d*P;
        Hxxp = blkdiag(Hxxp, H3d);
    end
    
    Hxx = Hxx + Hxxp;
      
end

function cost = cost_fct(X, U, param, R)
    %{
    X (NxM)

    %}
    X1 = reshape(X-param.Obj3d(:,1), (param.nbVarX)*param.nbData,1);
    Q1 = param.Q1;
    Q2 = param.Q2;
    cost = X1'*Q1*X1 + U'*R*U;
    
    % projecting points
    Xp = Point_Projection(X, param.ThetaX, param.ThetaY, param.ThetaZ, param.pov);
    X1p = reshape(Xp-param.Obj(:,1), (param.nbVarX-1)*param.nbData, 1);
    cost = cost + X1p'*param.Qp*X1p;
    
    for j = 1:param.nbObj-1             % iterating over points of the trajectory
       for i = 2:param.nbData-1         % iterating over the number of objectives
           
            X2 = reshape(Xp-param.traj(:,i,j), (param.nbVarX-1)*param.nbData,1);
            cost = cost - (X2'*Q2*X2);
            
        end
    end
    
end


function cost = cost_plot(X, param)

    Q2 = param.Q2(1:2,1:2);
    [~,M] = size(X);
    cost = zeros(1,M);
    
    for k = 1:M       
       for i = 1:param.nbData
           for j = 1:param.nbObj-1
               X2 = X(:,k) - param.traj(:,i,j);
               cost(:,k) = cost(:,k)-X2'*Q2*X2;
           end
       end
           
    end

end
