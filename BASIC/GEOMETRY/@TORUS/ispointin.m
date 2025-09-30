function [rep,P] = ispointin(C,P)
% function [rep,P] = ispointin(C,P)

c = getcoord(P);

% Radius and height
r = C.r;
h = C.h;

% Rotation matrix
v = [C.vx, C.vy];
n = [C.nx, C.ny, C.nz];
R = calcrotation(C,v,n);

% Apply inverse transform to point: project into the cylinder local frame
center = [C.cx, C.cy, C.cz];
vec = c - center; % vector from center to point
vec = vec * R';

% Height condition: point lies between base and top along cylinder axis
tol = getfemobjectoptions('tolerancepoint');
inHeight = (vec(:,3) >= -tol) & (vec(:,3) <= h + tol);

% Radial condition: point lies within cylinder radius
% d = sqrt(vec(:,1).^2 + vec(:,2).^2);
d = hypot(vec(:,1), vec(:,2));
inRadius = d <= r + eps;

% Both conditions: point lies between base and top along cylinder axis and within circle radius
rep = find(inHeight & inRadius);

if nargout==2
    P = P(rep);
end
