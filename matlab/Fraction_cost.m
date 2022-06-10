
function Fraction_cost

%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param.dt = 1E-2;        % Time step length
param.nbData = 50;      % Number of datapoints
param.nbIter = 5000;    % Number of iterations for iLQR
param.nbVarX = 2;       % State space dimension (x1,x2,x3)
param.nbVarU = 2;       % Control space dimension (dx1,dx2,dx3)
param.nbVarF = 2;       % Objective function dimension (f1,f2,f3, with f3 as orientation)
param.r = 0.1;          % Control weight term

% defining the two objectives
param.Obj = [[20.; -1.], [20.; 1.]];
%param.Obj = [[22.; -4.], [5.; 5.]];
%param.Obj = [[15.; 6.], [-2.; 3.]];
%param.Obj = [[23; 4], [24; 7]];

param.traj = zeros(param.nbVarX, param.nbData);

% generate the trajectory to be avoided
param.traj(1,:) = linspace(0,param.Obj(1,2),param.nbData);
param.traj(2,:) = linspace(0,param.Obj(2,2),param.nbData);

% computing the covariance matrix of the various goals
X_obj = param.Obj - mean(param.Obj, 2);
M = 2;
Cov = X_obj*X_obj'/(M-1);
epsilon = 1E-4 * diag(ones(M));
Cov = Cov + epsilon; % ensure not to have a singular matrixx

% Compute the eigenvectos of the covariance matrix
[eig_vec, D] = eig(Cov);
eig_max = max(D, [], 'all');

%tuning the weights as a function of the eigenvalues
%start by sorting the eigenvalues and rearanging eigenvectors
[~,I] = sort(diag(D),'descend');
eig_vec = eig_vec(:,I);
D = diag([0.1,1]);


scaling_fact = 1;

param.w2 = 100*ones(1,param.nbData);

q = zeros(1,param.nbData);
q(:) = scaling_fact;
param.Q2 = kron(diag(q), eig_vec*D*eig_vec');

% Q2 brings end effector to right goal
param.Qf = 10*param.w2(1)*param.nbData;
Mat1 = eye(param.nbData);
Mat1(param.nbData, param.nbData) = param.Qf; %terminal weight
param.Q1 = kron(Mat1, eig_vec*eye(param.nbVarX)*eig_vec');


% Defining the control matrix
R = speye((param.nbData-1)*param.nbVarU)*param.r;


%% Iterative LQR (iLQR)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u = zeros(param.nbVarU*(param.nbData-1), 1); %Initial commands
x0 = [0;0]; %Initial robot pose

%Transfer matrices (for linear system as single integrator)
Su0 = [zeros(param.nbVarX, param.nbVarX*(param.nbData-1)); tril(kron(ones(param.nbData-1), eye(param.nbVarX)*param.dt))];
Sx0 = kron(ones(param.nbData,1), eye(param.nbVarX));

for n=1:param.nbIter	
	x = reshape(Su0 * u + Sx0 * x0, param.nbVarX, param.nbData); %System evolution  

    %computing the Jacobians and Hessian of the cost function around the
    %current trajectory
    Jx = Gradient_x(x, param);
    Ju = 2*R*u;
    
    Hxx = Hess_x(x, param);
    Huu = 2*R;
    
    %calculating the update step
    du = (Su0'*Hxx*Su0 + Huu)^(-1)*(-Su0'*Jx-Ju);
    
	%Estimate step size with backtracking line search method
	
	alpha = 1;
	cost0 = cost_fct(x, u, param, R);
    
	while 1
		utmp = u + du * alpha;
		xtmp = reshape(Su0 * utmp + Sx0 * x0, param.nbVarX, param.nbData); %System evolution
        
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
if abs(x(1,param.nbData) - param.Obj(1,1)) < 1E-5
    if abs(x(2, param.nbData) -param.Obj(2,1)) < 1E-5
        
        display("reached goal point");
    end
end

    
figure('position',[10,10,550,550], 'color', [1 1 1]);
hold on;
axis off
%plotting the isolines of the cost function
Xs = linspace(-15,30,300);
[Grid_X, Grid_Y] = meshgrid(Xs);
Grid_X = reshape(Grid_X, 1, []);
Grid_Y = reshape(Grid_Y, 1, []);
Grid_coord = [Grid_X; Grid_Y];
costP = cost_plot(Grid_coord, param);
costP = reshape(costP, 300, 300);

size_rec = 0.5;
rec1_x = [param.Obj(1,1)-size_rec, param.Obj(1,1)-size_rec,  param.Obj(1,1)+size_rec, param.Obj(1,1)+size_rec];
rec1_y = [param.Obj(2,1)-size_rec, param.Obj(2,1)+size_rec, param.Obj(2,1)+size_rec, param.Obj(2,1)-size_rec];
rec2_x = [param.Obj(1,2)-size_rec, param.Obj(1,2)-size_rec,  param.Obj(1,2)+size_rec, param.Obj(1,2)+size_rec];
rec2_y = [param.Obj(2,2)-size_rec, param.Obj(2,2)+size_rec, param.Obj(2,2)+size_rec, param.Obj(2,2)-size_rec];


patch(rec1_x, rec1_y, [0.4660 0.6740 0.1880]);
patch(rec2_x, rec2_y, [0.8500 0.3250 0.0980]);

plot(x(1,:),x(2,:), "LineWidth", 1.5, "Marker", 'o');
xlim([-1,22])
ylim([-12.5,12.5])

%plotting fastest path
plot([x0(1) param.Obj(1,2)], [x0(2) param.Obj(2,2)], 'color', [0.4940 0.1840 0.5560], "LineWidth", 0.8, "LineStyle", '--');
plot([x0(1) param.Obj(1,1)], [x0(2) param.Obj(2,1)], 'color', [0.4940 0.1840 0.5560], "LineWidth", 0.8, "LineStyle", '--');


legend('Location', 'best')
legend("desired target", "undesired target", "legible trajectory", "optimal trajectory");


hold off;
figure
hold on
grid on
Xs = linspace(-15,30,300);
[Grid_X, Grid_Y] = meshgrid(Xs);
mesh(Grid_X, Grid_Y,costP)
patch(rec1_x, rec1_y, [0.4660 0.6740 0.1880]);
patch(rec2_x, rec2_y, [0.8500 0.3250 0.0980]);
view(3)
title("fractionary cost")
hold off


end



function Jx = Gradient_x(X, param)
    %{ 
    Computes the gradient of the cost function with respect to the x
    X: (NxM) datapoints
    Jx: (NMx1) vector containing the gradient of each step

    %}
    X1 = reshape(X-param.Obj(:,1), param.nbVarX*param.nbData,1);
  
    Jx = 2*param.Q1*X1;
    
    for k = 1:param.nbData
        Q2 = param.Q2(2*k-1:2*k, 2*k-1:2*k);
        for i = 2:param.nbData
            X2 = X(:,k)-param.traj(:,i);

            Jx(2*k-1:2*k) = Jx(2*k-1:2*k) - param.w2(i)*2*Q2*X2/((X2'*Q2*X2+1)^2);
        end
    end

end

function Jx = Gradient_plot(X, param)

    Q1 = param.Q1(1:2,1:2);
    Q2 = param.Q2(1:2,1:2);
    
    [~,M] = size(X);
    
    Jx = zeros(2,M);
    for k = 1:M
        Jx(:,k) = 2*Q1*(X(:,k)-param.Obj(:,1));
        
        for i = 2:param.nbData
            X2 = X(:,k)-param.traj(:,i);
        
            Jx(:,k) = Jx(:,k) - param.w2*2*Q2*X2/((X2'*Q2*X2)^2);
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
    
    Q1 = param.Q1;
   
    Hxx = 2*Q1;
    for k = 1:param.nbData
        Q2 = param.Q2(2*k-1:2*k,2*k-1:2*k);
        for i = 2:param.nbData
            X2 = X(:,k)-param.traj(:,i);

            Hxx(2*k-1:2*k,2*k-1:2*k) = Hxx(2*k-1:2*k,2*k-1:2*k)-param.w2(i)*((2*Q2/((X2'*Q2*X2+1)^2) - 8*Q2*X2*(Q2*X2)'/((X2'*Q2*X2+1)^3)));
        end
    end
     
    
end

function cost = cost_fct(X, U, param, R)
    %{
    X (NxM)

    %}
    X1 = reshape(X-param.Obj(:,1), param.nbVarX*param.nbData,1);
    Q1 = param.Q1;
   
    cost = X1'*Q1*X1 + U'*R*U;
    
    for k = 1:param.nbData
        Q2 = param.Q2(2*k-1:2*k,2*k-1:2*k);
        for i = 2:param.nbData
            

            X2 = X(:,k)-param.traj(:,i);
            cost = cost + param.w2(i)/(X2'*Q2*X2+1);

        end
    end
    
end

function cost = cost_plot(X, param)
    Q1 = param.Q1(1:2,1:2);
    Q2 = param.Q2(1:2,1:2);
    
    [~,M] = size(X);
    cost = zeros(1,M);
    for k = 1:M
       
        for i = 2:param.nbData
   
        
            X2 = X(:,k)-param.traj(:,i);
            cost(k) = cost(k) + param.w2(i)*0.01/(X2'*Q2*X2+1);

        end
        
    end
    


end
