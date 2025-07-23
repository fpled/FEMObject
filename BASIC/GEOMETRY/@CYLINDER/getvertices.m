function P = getvertices(C)
% function P = getvertices(C)

tol = getfemobjectoptions('tolerancepoint');

% Radius, height and opening angle
r = C.r;
h = C.h;
angle = C.angle;

% Vertices for full cylinder
% P = cell(1,8);
% P{1} = [ r,  0, 0]; % +x, base
% P{2} = [ 0,  r, 0]; % +y, base
% P{3} = [-r,  0, 0]; % -x, base
% P{4} = [ 0, -r, 0]; % -y, base
% P{5} = [ r,  0, h]; % +x, top
% P{6} = [ 0,  r, h]; % +y, top
% P{7} = [-r,  0, h]; % -x, top
% P{8} = [ 0, -r, h]; % -y, top

% Canonical order: (+x, +y, -x, -y)
base_angles = [0, pi/2, pi, 3*pi/2];
angles = base_angles(base_angles < angle + tol); % include only up to angle
if abs(angle - 2*pi) >= tol && angle > 0
    angles = [angles, angle];
end
% angles = 0; % always include start point at +x
% if angle > pi/2 - tol  , angles(end+1) = pi/2;   end
% if angle > pi - tol    , angles(end+1) = pi;     end
% if angle > 3*pi/2 - tol, angles(end+1) = 3*pi/2; end
% if abs(angle - 2*pi) >= tol && angle > 0
%     angles(end+1) = angle; % add endpoint at angle if not a full circle
% end

% Vertices at base and top
P = cell(1, 2*length(angles));
for i=1:length(angles)
    a = angles(i); % angle
    P{i} = [r*cos(a), r*sin(a), 0]; % base vertex
    P{i+length(angles)} = [r*cos(a), r*sin(a), h]; % top vertex
end

%% Old version
% Rotate around axis n = [nx, ny, nz] by angle of rotation phi =
% atan2(vy, vx) using tangent vector v = [vx, vy] via Rodrigues'
% rotation formula
%% New version
% Twist the XY plane about z = [0, 0, 1] by phi = atan2(vy, vx)
% using tangent vector v = [vx, vy], then tilt from z axis to normal
% vector n = [nx, ny, nz] so that the circle's normal is n regardless
% of v = [vx, vy]
v = [C.vx, C.vy];
n = [C.nx, C.ny, C.nz];
R = calcrotation(C,v,n);

% Translate to center c = [cx, cy, cz]
c = [C.cx, C.cy, C.cz];

for i=1:length(P)
    P{i} = P{i}*R + c;
end
