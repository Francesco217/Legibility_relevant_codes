function Projection = Point_projection_inverse(X, ThetaX, ThetaY, ThetaZ, pos, param, x0, fill)
%{
  X: (dimX x NbPoints) matrix containing all the points to be projected
  ThetaX/Y/Z: rotations correspoding to the angle of view of the obserever,
              The angle of view should be the result of rotation in the
              order ThetaX, Thetay and then ThetaZ

  pos: (dimX x 1) position of the observer

  return (dimX - 1 x NbPoints): projection of the NbPoints on the desired
                                plane we consider the last two coordinates
                                as the plane that intersts us

%}


Rx = [1 0 0; 0 cos(ThetaX) -sin(ThetaX); 0 sin(ThetaX) cos(ThetaX)];
Ry = [cos(ThetaY) 0 sin(ThetaY); 0 1 0; -sin(ThetaY) 0 cos(ThetaY)];
Rz = [cos(ThetaZ) -sin(ThetaZ) 0; sin(ThetaZ) cos(ThetaZ) 0; 0 0 1];
Rtot = Rz*Ry*Rx;
Rtot = inv(Rtot);
x0p = Rtot*x0 - Rtot*pos;
goal = Rtot*param.Obj3d(:,1) - Rtot*pos;

%filling the missing direction with the trajectory of a standard LQR
%controller
if strcmp(param.fill, 'lqr')
    A = 1;
    B = 1;
    Q = param.Q1(1,1)/200;
    R = param.r;
    [K,~,~] = dlqr(A,B,Q,R,0);

    x_fill = zeros(1,param.nbData);
    x_fill(1) = x0p(2);
    for i = 2:param.nbData
        x_fill(i) = (1-K)*(x_fill(i-1)-goal(2))+goal(2);  
        
    end
end

if strcmp(param.fill, 'linear')
    x_fill = linspace(x0p(2), goal(2), param.nbData);
end

X_ext = [X(1,:); x_fill ;X(2,:)];

Projection = X_ext + Rtot*pos;
Projection = Rtot'*Projection; 






end