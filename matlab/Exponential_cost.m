function Exponential_cost

%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param.dt = 1E-2;        %Time step length
param.nbData = 50;      %Number of datapoints
param.nbIter = 1500;    %Number of iterations for iLQR
param.nbVarX = 2;       %State space dimension (x1,x2,x3)
param.nbVarU = 2;       %Control space dimension (dx1,dx2,dx3)
param.nbVarF = 2;       %Objective function dimension (f1,f2,f3, with f3 as orientation)
param.r = 1;            %Control weight term
param.q2 = 1000;

% defining the two objectives
param.Obj = [[20; 1], [20; -1]];
%param.Obj = [[25;4], [23; 9]];
%param.Obj = [[20; 10], [-3; -7]];
%param.Obj = [[1;5], [-1;5]];

% computing the covariance matrix of the various goals
X_obj = param.Obj - mean(param.Obj, 2);
M = 2;
Cov = X_obj*X_obj'/(M-1);


% Compute the eigenvectos of the covariance matrix
[eig_vec, D] = eig(Cov);
eig_max = max(D, [], 'all');

% start by sorting the eigenvalues and rearanging eigenvectors
[~,I] = sort(D,2,'descend');
eig_vec = eig_vec(:,I(1,:));
Diag = diag([200,1]);
Cov = eig_vec*Diag*eig_vec';

epsilon = 1E-1 * diag(ones(M));
Cov = Cov + epsilon;% ensure not to have a singular matrix

% gererating the seed points for the avoidance of the trajectory
traj = zeros(param.nbVarX, param.nbData);
traj(1,:) = linspace(0,param.Obj(1,2), param.nbData);
traj(2,:) = linspace(0,param.Obj(2,2), param.nbData);

% Defining the Precision matrices 
% precision matrix Q1 is related to the cost of being close to the wrong goa
q = zeros(1,param.nbData);
q(1:10) = linspace(20,1,10);


%Q2 brings end effector to right goal
Mat1 = eye(param.nbData);
Mat1(param.nbData, param.nbData) = param.nbData*param.q2*10000; %terminal weight
Q1 = kron(Mat1, eye(param.nbVarX)*100);

%Defining the control matrix
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
    X1 = reshape(x-param.Obj(:,1), param.nbVarX*param.nbData,1);
    Jx = 2*Q1*X1 + dGauss(Cov, x,traj,param);
    Ju = 2*R*u;
    
    Hxx = 2*Q1 + Hess_Gauss(Cov, x, traj,param);
    Huu = 2*R;
    %calculating the update step
    du = (Su0'*Hxx*Su0 + Huu)^(-1)*(-Su0'*Jx-Ju);
    
 
	%Estimate step size with backtracking line search method
	
	alpha = 1;
	cost0 =  cost(x, traj, param, Cov, u, Q1,R);
   
	while 1
		utmp = u + du * alpha;
		xtmp = reshape(Su0 * utmp + Sx0 * x0, param.nbVarX, param.nbData); %System evolution 
		
        cost1 = cost(xtmp, traj, param, Cov, utmp, Q1,R);
        
		if cost1 < cost0 || alpha < 1E-3
			break;
		end
		alpha = alpha * 0.5;
	end
	u = u + du * alpha;
	
	if norm(du * alpha) < 1E-2
		break; %Stop iLQR when solution is reached
    end
    clc
    display(n)
end
disp(['iLQR converged in ' num2str(n) ' iterations.']);


%% Plot state space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('position',[10,10,700,700], 'color', [1 1 1]);
hold on;
axis off
size_rec = 0.75;
rec1_x = [param.Obj(1,1)-size_rec, param.Obj(1,1)-size_rec,  param.Obj(1,1)+size_rec, param.Obj(1,1)+size_rec];
rec1_y = [param.Obj(2,1)-size_rec, param.Obj(2,1)+size_rec, param.Obj(2,1)+size_rec, param.Obj(2,1)-size_rec];
rec2_x = [param.Obj(1,2)-size_rec, param.Obj(1,2)-size_rec,  param.Obj(1,2)+size_rec, param.Obj(1,2)+size_rec];
rec2_y = [param.Obj(2,2)-size_rec, param.Obj(2,2)+size_rec, param.Obj(2,2)+size_rec, param.Obj(2,2)-size_rec];
patch(rec1_x, rec1_y, [0.4660 0.6740 0.1880]);
patch(rec2_x, rec2_y, [0.8500 0.3250 0.0980]);

plot(x(1,:),x(2,:), "LineWidth", 1.5, "Marker", 'o')
xlim([-1,22])
ylim([-12.5,12.5])

%plotting fastest path
plot([x0(1) param.Obj(1,2)], [x0(2) param.Obj(2,2)], 'color', [0.4940 0.1840 0.5560], "LineWidth", 0.8, "LineStyle", '--');
plot([x0(1) param.Obj(1,1)], [x0(2) param.Obj(2,1)], 'color', [0.4940 0.1840 0.5560], "LineWidth", 0.8, "LineStyle", '--');
legend('Location', 'best')
legend("desired target", "undesired target", "legible trajectory", "optimal trajectory");


scatter(traj(1,:),traj(2,:))

hold off;

figure
Xs = linspace(-10,30,300);
[Grid_X, Grid_Y] = meshgrid(Xs);
Grid_X = reshape(Grid_X, 1, []);
Grid_Y = reshape(Grid_Y, 1, []);
Grid_coord = [Grid_X; Grid_Y];
costP = cost_plot(Grid_coord, param, traj, Cov);
costP = reshape(costP, 300, 300);
[Grid_X, Grid_Y] = meshgrid(Xs);
mesh(Grid_X, Grid_Y,costP)



end



function G = Gaussian(Sigma, x,traj,param)
    %{
       sigma: covariance matrix
       x : (dimX x nbVar) 
       traj : (dimX x nbVar) trajectory points
       
        G is (1xM) row vector

    %}

    %reshape the data in a convenient format
    X_diff = x-traj;
    
    %creating the resulting value of the exponential part of the Gauss function
    exp_arg = diag(exp(-0.5*X_diff'*inv(Sigma)*X_diff));
    G = exp_arg;

end

function dPDF = dGauss(Sigma, X,traj,param)
%{
    X is (NxM) where M is the number of points and N is the dimension
    dPDF is the derivative of the Gaussian distribution computed on each of
    the points and is given as a M*Nx1 long vector
%}
    dPDF = zeros(param.nbVarX*param.nbData,1);
    Sigma_inv = inv(Sigma);
    
    for k = 1:param.nbData      %iterate over the points of optimisation
        for i = 2:param.nbData  %iterate over the trajectory to avoid
            
            X_diff = X(:,k)-traj(:,i);
            dPDF(2*k-1:2*k,:) = dPDF(2*k-1:2*k,:) -param.q2*Sigma_inv*(X_diff)*Gaussian(Sigma, X_diff, traj(:,i), param);
        end
    end
    %reshaping it as a clomumn vector
   % dPDF = dPDF';
end

function H_gauss = Hess_Gauss(Sigma, X, traj,param)
    %{
        X is (NxM), the M datapoints
        Sigma (NxN) covariance matrix
        traj (NxM), trajectory to avoid
        param structure of parameter
        
        H is (NMxNM) the Hessian whose diag is the hessian for each
        datapoint
    %}

    
    Sigma_inv = inv(Sigma);
    H_gauss = zeros(param.nbVarX*param.nbData);
    
    for k = 1:param.nbData
        for i = 2:param.nbData
            %reshaping the data in a long vector
            X_diff = X(:,k)-traj(:,i);

            %computing the product X_diff*X_diff and taking only the diagonal
            X_diag = X_diff*X_diff'; %extracing the diagonal
           

            %creating the exponential term of the Hessian matrix in the correct
            %format
            exp_term = Gaussian(Sigma, X_diff, traj(:,i), param);

            %computing the Hessian matrix
            H_gauss(2*k-1:2*k,2*k-1:2*k) = H_gauss(2*k-1:2*k,2*k-1:2*k) - param.q2*Sigma_inv*(1-X_diag*Sigma_inv)*exp_term;
            
        end
    end
    


   
end

function c = cost(X, traj, param, Sigma, U, Q,R)
    
    %reshape into a long 1D vector
    X_diff1 = reshape(X-param.Obj(:,1), param.nbVarX*param.nbData,1);
    
    %attractive term of the cost function
    
    c = X_diff1'*Q*X_diff1 + U'*R*U;
    for i = 1:param.nbData
        
        c = c  + param.q2*sum(Gaussian(Sigma, X, traj(:,i), param));
    end
    
end

function cost = cost_plot(X, param, traj, Sigma)
  
    
    [~,M] = size(X);
    cost = zeros(1,M);
    
    for k = 1:M
        for i = 1:param.nbData
            cost(:,k) = cost(:,k) + param.q2*Gaussian(Sigma, X(:,k), traj(:,i), param);
        end
    end
  

end
