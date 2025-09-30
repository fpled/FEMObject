function cminor = getcminor(T)
% function cminor = getcminor(T)

tol = getfemobjectoptions('tolerancepoint');

% Major radius and opening angle
r1 = T.r1;
angle = T.angle;
if isstring(angle), angle = char(angle); end
if ischar(angle),   angle = str2num(lower(angle)); end

isfull = abs(angle - 2*pi) < tol;

if isfull
    % Full torus
    cminor = [ r1,   0, 0;  % +x
                0,  r1, 0;  % +y
              -r1,   0, 0;  % -x
                0, -r1, 0]; % -y
else
    % Partial torus
    base_angles = [0, pi/2, pi, 3*pi/2]; % canonical angles
    angles = base_angles(base_angles + tol < angle); % keep only canonical angles strictly before angle (within tol) 
    if angle > tol && (isempty(angles) || abs(angle - angles(end)) > tol)
        angles(end+1) = angle; % append endpoint at angle
    end
    N = numel(angles);
    
    % Vertices at each major angle
    cminor = zeros(N,3);
    for i=1:N
        t = angles(i);
        [ct,st] = cos_sin_snap(t,tol); % snap to canonical to avoid roundoff
        cminor(i,:) = [r1*ct, r1*st, 0];
    end
end

% Rotation matrix
v = [T.vx, T.vy];
n = [T.nx, T.ny, T.nz];
R = calcrotation(T,v,n);

% Center
c = [T.cx, T.cy, T.cz];

% Rotate into global frame and translate to center
cminor = cminor*R + c;

end

function [ct,st] = cos_sin_snap(t,tol)
if abs(t-0) < tol,        ct = 1;      st = 0;  % +x
elseif abs(t-pi/2) < tol, ct = 0;      st = 1;  % +y
elseif abs(t-pi) < tol,   ct = -1;     st = 0;  % -x
elseif abs(t-3*pi/2)<tol, ct = 0;      st = -1; % -y
else,                     ct = cos(t); st = sin(t);
end
end
