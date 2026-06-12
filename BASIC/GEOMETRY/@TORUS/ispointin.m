function [rep,P] = ispointin(T,P)
% function [rep,P] = ispointin(T,P)

tol = getfemobjectoptions('tolerancepoint');

c = getcoord(P);

% Major radius, minor radius, and opening angle
r1 = T.r1;
r2 = T.r2;
angle = T.angle;

if isstring(angle), angle = char(angle); end
if ischar(angle),   angle = str2num(lower(angle)); end

isfull = abs(angle - 2*pi) < tol;

% Rotation matrix from torus local frame to global frame
v = [T.vx, T.vy];
n = [T.nx, T.ny, T.nz];
R = calcrotation(T,v,n);

% Apply inverse transform to point: project into torus local frame
center = [T.cx, T.cy, T.cz];
vec = c - center; % vector from center to point
vec = vec * R';

% Torus local coordinates
% d = sqrt(vec(:,1).^2 + vec(:,2).^2); % radial distance
d = hypot(vec(:,1),vec(:,2));        % radial distance
z = vec(:,3);                        % axial coordinate

% Tube condition: point lies within tube
inTube = (d - r1).^2 + z.^2 <= (r2 + tol)^2;

% Angular condition for partial torus
if isfull
    inAngle = true(size(inTube));
else
    theta = atan2(vec(:,2),vec(:,1));
    theta(theta < 0) = theta(theta < 0) + 2*pi;

    tolAngle = tol / max(r1,eps);
    inAngle = (theta <= angle + tolAngle) | ...
              (theta >= 2*pi - tolAngle);
end

% Both conditions: point lies within tube with angular opening (only for partial torus)
rep = find(inTube & inAngle);

if nargout==2
    P = P(rep);
end
