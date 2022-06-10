
function Projection = Point_Projection(X, ThetaX, ThetaY, ThetaZ, pos)
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

Rtot = inv(Rz*Ry*Rx); % take inverse because object turns other direction than direction camera turns

Projection = Rtot*X - Rtot*pos;

Projection = Projection([1,3],:);


end