function P = getvertices(T)
% function P = getvertices(T)

tol = getfemobjectoptions('tolerancepoint');

% Major radius, minor radius, and opening angle
r1 = T.r1;
r2 = T.r2;
angle = T.angle;
if isstring(angle), angle = char(angle); end
if ischar(angle),   angle = str2num(lower(angle)); end

isfull = abs(angle - 2*pi) < tol;

if isfull
    % Full torus
    % Vertices per major angle order: outer, top, inner, bottom
    P = cell(1,16);
    c1 = [r1, 0, 0]; % +x (t = 0)
    P{1} = c1 + [ r2, 0,   0]; % outer
    P{2} = c1 + [  0, 0,  r2]; % top
    P{3} = c1 + [-r2, 0,   0]; % inner
    P{4} = c1 + [  0, 0, -r2]; % bottom
    
    c2 = [0, r1, 0]; % +y (t = pi/2)
    P{5} = c2 + [0,  r2,   0]; % outer
    P{6} = c2 + [0,   0,  r2]; % top
    P{7} = c2 + [0, -r2,   0]; % inner
    P{8} = c2 + [0,   0, -r2]; % bottom
    
    c3 = [-r1, 0, 0]; % -x (t = pi)
    P{9}  = c3 + [-r2, 0,   0]; % outer
    P{10} = c3 + [  0, 0,  r2]; % top
    P{11} = c3 + [ r2, 0,   0]; % inner
    P{12} = c3 + [  0, 0, -r2]; % bottom
    
    c4 = [0, -r1, 0]; % -y (t = 3*pi/2)
    P{13} = c4 + [0, -r2,   0]; % outer
    P{14} = c4 + [0,   0,  r2]; % top
    P{15} = c4 + [0,  r2,   0]; % inner
    P{16} = c4 + [0,   0, -r2]; % bottom
else
    % Partial torus
    base_angles = [0, pi/2, pi, 3*pi/2]; % canonical angles
    angles = base_angles(base_angles + tol < angle); % keep only canonical angles strictly before angle (within tol)
    if angle > tol && (isempty(angles) || abs(angle - angles(end)) > tol)
        angles(end+1) = angle; % append endpoint at angle
    end
    N = numel(angles);
    
    % Vertices per angle order: outer, top, inner, bottom
    P = cell(1,4*N);
    for i=1:N
        t = angles(i);
        [ct,st] = cos_sin_snap(t,tol); % snap to canonical to avoid roundoff
        k = 4*(i-1);
        P{k+1} = [(r1+r2)*ct, (r1+r2)*st,   0]; % outer
        P{k+2} = [     r1*ct,      r1*st,  r2]; % top
        P{k+3} = [(r1-r2)*ct, (r1-r2)*st,   0]; % inner
        P{k+4} = [     r1*ct,      r1*st, -r2]; % bottom
    end
end

% Rotation matrix
v = [T.vx, T.vy];
n = [T.nx, T.ny, T.nz];
R = calcrotation(T,v,n);

% Center
c = [T.cx, T.cy, T.cz];

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