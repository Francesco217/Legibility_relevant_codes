function q = end_position_with_proj(param)
%{
    This function computes the target end position in configuration space
    
    param: structure containing the problem's parameters

    q_N: (param.nbVarX x 1) computed goal end position.

%}

nbIter = 400;

% attractive matrix
param2.Q1 = zeros(param.nbVarF, param.nbVarF, length(param.idx));

for k = 1:length(param.idx)
   if k == length(param.idx)
        param2.Q1(:,:,k) = eye(param.nbVarF)*10000;
   else
       param2.Q1(:,:,k) = eye(param.nbVarF)*10;
   end
end

% computing the covariance matrix of the various goals in projection space
X_obj = param.Obj - mean(param.Obj, 2);

Cov = X_obj*X_obj'/(param.nbObj-1);
epsilon = 1E-4 * diag(ones(1,param.nbVarF-1));
Cov = Cov+epsilon; % ensure not to have a singular matrixx

% extract the eigenvectors 
% Compute the eigenvectos of the covariance matrix
[eig_vec, D] = eig(Cov);
eig_max = max(D, [], 'all');

% tuning the weights as a function of the eigenvalues
% start by sorting the eigenvalues and rearanging eigenvectors
[~,I] = sort(diag(D),'descend');
eig_vec = eig_vec(:,I);

D = diag([10,0.1]);

param2.Q2 = eig_vec*D*eig_vec'*(10);

% finding the soltion of the problem using gradient descent
q = zeros(param.nbVarX,1);

for e = 1:nbIter
    
    gamma = 0.5;
    c0 = cost_fct(q, param, param2);
    d = grad(q, param, param2);
    dq = d;
    % compute the update step
    qtmp = q;
    while cost_fct(qtmp, param, param2) > c0 || gamma > 1E-3
        qtmp = q - gamma*dq;
        gamma = 0.5*gamma;
       
    end
    q = qtmp;
    if abs(cost_fct(q, param, param2) - c0) < 1E-4
        break;
    end
end

display(['end position computed']) 

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
	J = logmap(F2, F1) / e; % Error by considering manifold
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
	e0 = [0; 0; 0; 1]; %Origin of the manifold
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




e = 1E-4;
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

function  cost = cost_fct(q, param, param2)
    
   % start by computing the forward kinematics of all the articulations
   [F,~] = fkin0(q, param);
   F = F(1:3,param.idx);
   Fp = Point_Projection(F, param.ThetaX, param.ThetaY, param.ThetaZ, param.pov);
   cost = 0;
   for k = 1:length(param.idx)
      X1 = F(:,k) - param.Obj3d(:,1);
      
      cost = cost + X1'*param2.Q1(:,:,k)*X1;
      
      if not(k == length(param.idx))
        for n = 2:param.nbObj
            X2 = Fp(:,k) - param.Obj(:,n);
            cost = cost - X2'*param2.Q2*X2;
        end
      end
       
   end

end

% Gradient calcuation -----------------------------------------------------

function J = grad(q, param, param2)
%{
    computes the gradient of the cost function at the current point

%}

% constructing the projection matrix
Rx = [1 0 0; 0 cos(param.ThetaX) -sin(param.ThetaX); 0 sin(param.ThetaX) cos(param.ThetaX)];
Ry = [cos(param.ThetaY) 0 sin(param.ThetaY); 0 1 0; -sin(param.ThetaY) 0 cos(param.ThetaY)];
Rz = [cos(param.ThetaZ) -sin(param.ThetaZ) 0; sin(param.ThetaZ) cos(param.ThetaZ) 0; 0 0 1];
P = inv(Rz*Ry*Rx);
P = P([1,3],:);

[F,~] = fkin0(q, param);
F = F(1:3,param.idx);
Fp = Point_Projection(F, param.ThetaX, param.ThetaY, param.ThetaZ, param.pov);

Jkin = Jkin_all(q,param);
J= zeros(param.nbVarX,1);
for k = 1:length(param.idx)
   X1 = F(:,k) - param.Obj3d(:,1);
   J1 = 2*param2.Q1(:,:,k)*(X1);
   
   if not(k == length(param.idx))
       for n = 2:param.nbObj
          X2 = Fp(:,k) - param.Obj(:,n); 
          J1 = J1 - 2*P'*param2.Q2*X2;
       end
   end
   
   J = J + Jkin(:,:,k)'*J1;
    
end

end

