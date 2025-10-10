function [rep,P] = ispointin(E,P)
% function [rep,P] = ispointin(E,P)

tol = getfemobjectoptions('tolerancepoint');

coord = getcoord(P);

% Semi-axes
a = E.a;
b = E.b;
c = E.c;

% Rotation matrix
v = [E.vx, E.vy];
n = [E.nx, E.ny, E.nz];
R = calcrotation(E,v,n);

% Center
center = [E.cx, E.cy, E.cz];

% Apply inverse transform to point: project into the ellipsoid local frame
vec = coord - center; % vector from center to point
vec = vec * R';

% Ellipsoidal condition: point lies within ellipsoid semi-axes
xc = vec(:,1)/a;
yc = vec(:,2)/b;
zc = vec(:,2)/c;
inEllipsoid = (xc.^2 + yc.^2 + zc.^2) <= 1 + tol;
rep = find(inEllipsoid);

if nargout==2
    P = P(rep);
end
