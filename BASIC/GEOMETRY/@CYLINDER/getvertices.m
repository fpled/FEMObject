function P = getvertices(C)
% function P = getvertices(C)

tol = getfemobjectoptions('tolerancepoint');

% Radius, height and opening angle
r = C.r;
h = C.h;
angle = C.angle;
if isstring(angle), angle = char(angle); end
if ischar(angle),   angle = str2num(lower(angle)); end

isfull = abs(angle - 2*pi) < tol;

if isfull
    % Full cylinder
    % Vertices at base (z = 0) and top (z = h)
    P = cell(1,8);
    % base (z = 0)
    P{1} = [ r,  0, 0]; % +x
    P{2} = [ 0,  r, 0]; % +y
    P{3} = [-r,  0, 0]; % -x
    P{4} = [ 0, -r, 0]; % -y
    % top (z = h)
    P{5} = [ r,  0, h]; % +x
    P{6} = [ 0,  r, h]; % +y
    P{7} = [-r,  0, h]; % -x
    P{8} = [ 0, -r, h]; % -y
else
    % Partial cylinder
    base_angles = [0, pi/2, pi, 3*pi/2]; % canonical angles
    angles = base_angles(base_angles + tol < angle); % keep only canonical angles strictly before angle (within tol) 
    if angle > tol && (isempty(angles) || abs(angle - angles(end)) > tol)
        angles(end+1) = angle; % append endpoint at angle
    end
    N = numel(angles);
    
    % Vertices at base (z = 0) and top (z = h)
    P = cell(1,2*N);
    for i=1:N
        t = angles(i);
        [ct,st] = cos_sin_snap(t,tol); % snap to canonical to avoid roundoff
        x = r * ct;
        y = r * st;
        P{i}   = [x, y, 0]; % base (z = 0)
        P{N+i} = [x, y, h]; % top  (z = h)
    end
end

% Rotation matrix
v = [C.vx, C.vy];
n = [C.nx, C.ny, C.nz];
R = calcrotation(C,v,n);

% Center
c = [C.cx, C.cy, C.cz];

% Rotate into global frame and translate to center
for i=1:numel(P)
    P{i} = P{i}*R + c;
end

end

function [ct,st] = cos_sin_snap(t,tol)
if abs(t-0) < tol,        ct = 1;      st = 0;  % +x
elseif abs(t-pi/2) < tol, ct = 0;      st = 1;  % +y
elseif abs(t-pi) < tol,   ct = -1;     st = 0;  % -x
elseif abs(t-3*pi/2)<tol, ct = 0;      st = -1; % -y
else,                     ct = cos(t); st = sin(t);
end
end