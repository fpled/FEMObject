function [rep,P] = ispointin(C,P)
% function [rep,P] = ispointin(C,P)

tol = getfemobjectoptions('tolerancepoint');

c = getcoord(P);

% Radius and height, and opening angle
r = C.r;
h = C.h;
angle = C.angle;

if isstring(angle), angle = char(angle); end
if ischar(angle),   angle = str2num(lower(angle)); end

isfull = abs(angle - 2*pi) < tol;

% Rotation matrix from cylinder local frame to global frame
v = [C.vx, C.vy];
n = [C.nx, C.ny, C.nz];
R = calcrotation(C,v,n);

% Apply inverse transform to point: project into cylinder local frame
center = [C.cx, C.cy, C.cz];
vec = c - center; % vector from center to point
vec = vec * R';

% Cylinder local coordinates
% d = sqrt(vec(:,1).^2 + vec(:,2).^2); % radial distance
d = hypot(vec(:,1),vec(:,2));        % radial distance
z = vec(:,3);                        % axial coordinate from base

% Height condition: point lies between base and top along cylinder axis
inHeight = (z >= -tol) & (z <= h + tol);

% Radial condition: point lies within cylinder radius
inRadius = d <= r + tol;

% Angular condition for partial cylinder
if isfull
    inAngle = true(size(d));
else
    theta = atan2(vec(:,2),vec(:,1));
    theta(theta < 0) = theta(theta < 0) + 2*pi;

    tolAngle = tol / max(r,eps);
    inAngle = (theta <= angle + tolAngle) | ...
              (theta >= 2*pi - tolAngle);
end

% Final condition: point lies between base and top along cylinder axis and
% within circle radius with opening angle (only for partial cylinder)
rep = find(inHeight & inRadius & inAngle);

if nargout==2
    P = P(rep);
end
