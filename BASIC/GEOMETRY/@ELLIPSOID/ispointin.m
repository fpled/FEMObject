function [rep,P] = ispointin(E,P)
% function [rep,P] = ispointin(E,P)

coord = getcoord(P);

% Semi-axes
a = E.a;
b = E.b;
c = E.c;

%% Old version
% Rotate around axis n = [nx, ny, nz] by angle of rotation phi =
% atan2(vy, vx) using tangent vector v = [vx, vy] via Rodrigues'
% rotation formula
%% New version
% Twist the XY plane about z = [0, 0, 1] by phi = atan2(vy, vx)
% using tangent vector v = [vx, vy], then tilt from z axis to normal
% vector n = [nx, ny, nz] so that the circle's normal is n regardless
% of v = [vx, vy]
v = [E.vx, E.vy];
n = [E.nx, E.ny, E.nz];
R = calcrotation(E,v,n);

% Apply inverse transform to point: project into the ellipsoid local frame
center = [E.cx, E.cy, E.cz];
vec = coord - center; % vector from center to point
vec = vec * R';

% Ellipsoidal condition: point lies within ellipsoid semi-axes
xc = vec(:,1)/a;
yc = vec(:,2)/b;
zc = vec(:,2)/c;
inEllipsoid = (xc.^2 + yc.^2 + zc.^2) <= 1 + eps;
rep = find(inEllipsoid);

if nargout==2
    P = P(rep);
end
