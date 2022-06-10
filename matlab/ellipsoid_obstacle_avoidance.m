
function ellipsoid_obstacle_avoidance

% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param.dt = 1E-2;        %Time step length
param.nbData = 50;      %Number of datapoints
param.nbIter = 1000;    %Number of iterations for iLQR
param.nbVarX = 3;       %State space dimension (x1,x2,x3)
param.nbVarU = 3;       %Control space dimension (dx1,dx2,dx3)
param.nbVarF = 2;       %Objective function dimension (f1,f2,f3, with f3 as orientation)
param.r = 100;          %Control weight term
param.nbObj = 2;
param.d = 4;             %distance for radius of repulsive cost

% defining the two objectives
%param.Obj3d = [[20.; -1.;0], [20.; 1.;0], [20;0;-1]];
%param.Obj3d = [[22.; -4.], [5.; 5.]];
%param.Obj3d = [[15.; 6.], [-2.; 3.]];
param.Obj3d = [[24; 7;5],[23; 4;10]];

% Initial robot pose
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
grid on;
view(3);
hold off

% ------------------------------------------------------------------------

% defining the point of view of the observer
param.pov = [15; -10; 5];
param.ThetaX = 0;
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
D = diag([10,5])/eig_max;


% TUNING PARAMETERS ------------------------------------------------------

% Q2 matrix, keep away from unwanted trajectories
scaling_fact = 1000;

param.w2 = linspace(3000, 1500, param.nbData);

q = zeros(1,param.nbData);
q(:) = scaling_fact;
param.Q2 = kron(diag(q), eig_vec*D*eig_vec');


% Q1 brings end effector to right goal
param.Qf = param.nbData*param.w2(1)*scaling_fact*100;
Mat1 = eye(param.nbData);
Mat1(param.nbData, param.nbData) = param.Qf; %terminal weight
param.Q1 = kron(Mat1, eye(param.nbVarX)*100);

% Qp brings end effector to right goal in projection space
Mat2 = eye(param.nbData);
Mat2(param.nbData, param.nbData) = param.Qf;
param.Qp = kron(Mat2, eig_vec*eye(param.nbVarX-1)*eig_vec'*0);


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
    
    % computing the Jacobians and Hessian of the cost function 
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
		xtmp = reshape(Su0 * utmp + Sx0 * x0_3d, param.nbVarX, param.nbData); %System evolution
        
        cost = cost_fct(xtmp, utmp, param, R);
        
        if cost < cost0 || alpha < 1E-4
			break;
        end
        
		alpha = alpha * 0.5;
    end
    
	u = u + du * alpha;
    
    clc;
    display(n);
	
	if norm(du * alpha) < 1E-2
		break; %Stop iLQR when solution is reached
    end
    
end
disp(['iLQR converged in ' num2str(n) ' iterations.']);

% retransforming everything in original frame

if abs(x3d(1,param.nbData) - param.Obj(1,1)) < 1E-5
    if abs(x3d(2, param.nbData) -param.Obj(2,1)) < 1E-5
        if abs(x3d(2, param.nbData)-param.Obj(2,1)) < 1E-5
        
            display("reached goal point");
        end
    end
end

%% Plot state space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot projection ---------------------------------------------------------
x = Point_Projection(x3d, param.ThetaX, param.ThetaY, param.ThetaZ, param.pov);

figure('position',[10,10,600,600], 'color', [1 1 1]);
hold on;


size_rec = 0.75;
rec1_x = [param.Obj(1,1)-size_rec, param.Obj(1,1)-size_rec,  param.Obj(1,1)+size_rec, param.Obj(1,1)+size_rec];
rec1_y = [param.Obj(2,1)-size_rec, param.Obj(2,1)+size_rec, param.Obj(2,1)+size_rec, param.Obj(2,1)-size_rec];

patch(rec1_x, rec1_y, [0.4660 0.6740 0.1880]);
for j = 2:param.nbObj
    rec2_x = [param.Obj(1,j)-size_rec, param.Obj(1,j)-size_rec, param.Obj(1,j)+size_rec, param.Obj(1,j)+size_rec];
    rec2_y = [param.Obj(2,j)-size_rec, param.Obj(2,j)+size_rec, param.Obj(2,j)+size_rec, param.Obj(2,j)-size_rec];
    patch(rec2_x, rec2_y, [0.8500 0.3250 0.0980]);
end
plot(x(1,:),x(2,:),"LineWidth", 1.5, "Marker", 'o')

%plotting fastest path
plot([x0(1) param.Obj(1,2)], [x0(2) param.Obj(2,2)], 'color', [0.4940 0.1840 0.5560], "LineWidth", 0.8, "LineStyle", '--');
plot([x0(1) param.Obj(1,1)], [x0(2) param.Obj(2,1)], 'color', [0.4940 0.1840 0.5560], "LineWidth", 0.8, "LineStyle", '--');


% grid on
axis off
hold off;

% Plotting the 3D trajectory-----------------------------------------------

figure('position',[10,10,600,600], 'color', [1 1 1]);
hold on
grid on

plot3(x3d(1,:),x3d(2,:),x3d(3,:),"LineWidth", 1.5);

% plotting the objectives
draw_cube(param.Obj3d(:,1), size_rec, [0.4660 0.6740 0.1880]);

for i = 2:param.nbObj
    
   draw_cube(param.Obj3d(:,i), size_rec, [0.8500 0.3250 0.0980]); 
    
end

plot3([x0_3d(1) param.Obj3d(1,1)], [x0_3d(2) param.Obj3d(2,1)], [x0_3d(3) param.Obj3d(3,1)], 'color', [0.4940 0.1840 0.5560], "LineWidth", 0.8, "LineStyle", '--');
plot3([x0_3d(1) param.Obj3d(1,2)], [x0_3d(2) param.Obj3d(2,2)], [x0_3d(3) param.Obj3d(3,2)], 'color', [0.4940 0.1840 0.5560], "LineWidth", 0.8, "LineStyle", '--');

Rx = [1 0 0; 0 cos(param.ThetaX) -sin(param.ThetaX); 0 sin(param.ThetaX) cos(param.ThetaX)];
Ry = [cos(param.ThetaY) 0 sin(param.ThetaY); 0 1 0; -sin(param.ThetaY) 0 cos(param.ThetaY)];
Rz = [cos(param.ThetaZ) -sin(param.ThetaZ) 0; sin(param.ThetaZ) cos(param.ThetaZ) 0; 0 0 1];

ey = [0;6;0];
ey_new = Rz*Ry*Rx*ey;

quiver3(param.pov(1), param.pov(2), param.pov(3), ey_new(1), ey_new(2), ey_new(3),"LineWidth", 5, 'color', 'b');
draw_cube(param.pov, 1, [0.9290 0.6940 0.1250])

xlim([-2,28]);
ylim([-15,15]);
zlim([-2,28]);
xlabel('x');
ylabel('y');
zlabel('z');

view(3);
hold off;


end

%GRADIENT function --------------------------------------------------------

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
    Jx = 2*param.Q1*X1;
    
    % Projecting X on 2D viewplane to calcultate 2D Jacobian in viewplane
    Xp = Point_Projection(X, param.ThetaX, param.ThetaY, param.ThetaZ, param.pov);
    % Xp = reshape(Xp, (param.nbVarX-1)*param.nbData, 1);
    X1p = reshape(Xp-param.Obj(:,1), (param.nbVarX-1)*param.nbData, 1);
    
    Jxp = 2*param.Qp*X1p;
    
    for k = 1:param.nbData          % iterate over the points of trajectorx
        Q2 = param.Q2(2*k-1:2*k,2*k-1:2*k);
        
        for j = 1:param.nbObj-1     % iterate over the number of different objectives
            for i = 2:param.nbData  % iterate over the trajectory to avoid

                dif1 = abs(Xp(1,k)-param.traj(1,i,j));
                dif2 = abs(Xp(2,k)-param.traj(2,i,j));
                j1 = 0;
                j2 = 0;
                f1 = 0;
                f2 = 0;
                if dif1 < param.d
                    f1 = param.d - dif1;
                    if Xp(1,k) > param.traj(1,i,j)
                        j1 = -1;
                    else
                        j1 = 1;
                    end
                    
                end
                
                if dif2 < param.d
                    f2 = param.d - dif2;
                    if Xp(2,k) > param.traj(2,i,j)
                        j2 = -1;
                    else
                        j2 = 1;
                    end
                end
                
                f = [f1;f2];
                Jac = diag([j1,j2]);
                
                Jxp(2*k-1:2*k,:) = Jxp(2*k-1:2*k,:) + 2*Jac'*Q2*f;
                
            end
        end   
    end
    
    % going from 2d Jacobian to the 3D one
    Jxp = reshape(Jxp, param.nbVarX-1, param.nbData);
    Jxp = reshape(P'*Jxp, param.nbVarX*param.nbData,1);
    Jx = Jx + Jxp;
 

end

% APRROX OF HESSIAN CALCULATION TO USE GAUSS NEWTON UPDATE STEP

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
    Hxxp = []; %H essian for all points 
    
    for k = 1:param.nbData
        H2d = 2*param.Qp(2*k-1:2*k, 2*k-1:2*k);
        Q2 = param.Q2(2*k-1:2*k,2*k-1:2*k);
        
        for j = 1:param.nbObj-1         % iterate over the number of different objectives
            for i = 2:param.nbData      % iterate over the trajectory    
            
                dif1 = abs(Xp(1,k)-param.traj(1,i,j));
                dif2 = abs(Xp(2,k)-param.traj(2,i,j));
                j1 = 0;
                j2 = 0;
  
                if dif1 < param.d
                    
                    if Xp(1,k) > param.traj(1,i,j)
                        j1 = -1;
                    else
                        j1 = 1;
                    end
                    
                end
                
                if dif2 < param.d
                    
                    if Xp(2,k) > param.traj(2,i,j)
                        j2 = -1;
                    else
                        j2 = 1;
                    end
                end
                
               
                Jac = diag([j1,j2]);
                
     
                H2d = H2d + 2*Jac'*Q2*Jac;
            end
        end
        
        % from projection space Hessian compute original space
        % Hessian

        H3d = P'*H2d*P;
        Hxxp = blkdiag(Hxxp, H3d);
    end
  
   Hxx = Hxx + Hxxp;
      
end


% COST FUNCTION -----------------------------------------------------------

function cost = cost_fct(X, U, param, R)
    %{
    X (NxM)

    %}
    X1 = reshape(X-param.Obj3d(:,1), (param.nbVarX)*param.nbData,1);
    Q1 = param.Q1;
    
    % attractive part of the cost
    cost = X1'*Q1*X1 + U'*R*U;
    
    % projecting points
    Xp = Point_Projection(X, param.ThetaX, param.ThetaY, param.ThetaZ, param.pov);
    X1p = reshape(Xp-param.Obj(:,1), (param.nbVarX-1)*param.nbData, 1);
    cost = cost + X1p'*param.Qp*X1p;
    
    for k = 1:param.nbData
        Q2 = param.Q2(2*k-1:2*k,2*k-1:2*k);
        for j = 1:param.nbObj-1        % iterating over the targets
           for i = 2:param.nbData      % iterating over the point of trajectory
                dif1 = abs(Xp(1,k)-param.traj(1,i,j));
                dif2 = abs(Xp(2,k)-param.traj(2,i,j));
                f1 = 0;
                f2 = 0;
                if dif1 < param.d
                    f1 = param.d - dif1;
                end
                if dif2 < param.d
                    f2 = param.d - dif2;
                end
                
                f = [f1;f2];
                cost = cost + f'*Q2*f;
               
                
            end
        end
    end
    
end

