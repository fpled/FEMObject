function P = getvertices(C)
% function P = getvertices(C)

% Radius and height
r = C.r;
h = C.h;

P = cell(1,8);
% Vertices
P{1} = [-r,  0, 0];
P{2} = [ 0, -r, 0];
P{3} = [ r,  0, 0];
P{4} = [ 0,  r, 0];
P{5} = [-r,  0, h];
P{6} = [ 0, -r, h];
P{7} = [ r,  0, h];
P{8} = [ 0,  r, h];

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

for i=1:8
    P{i} = P{i}*R + c;
end
