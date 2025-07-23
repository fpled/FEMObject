function P = getvertices(T)
% function P = getvertices(T)

% Radii
r1 = T.r1;
r2 = T.r2;

P = cell(1,8);
% Vertices
P{1} = [ r,  0, 0]; % +x, base
P{2} = [ 0,  r, 0]; % +y, base
P{3} = [-r,  0, 0]; % -x, base
P{4} = [ 0, -r, 0]; % -y, base
P{5} = [ r,  0, h]; % +x, top
P{6} = [ 0,  r, h]; % +y, top
P{7} = [-r,  0, h]; % -x, top
P{8} = [ 0, -r, h]; % -y, top

%% Old version
% Rotate around axis n = [nx, ny, nz] by angle of rotation phi =
% atan2(vy, vx) using tangent vector v = [vx, vy] via Rodrigues'
% rotation formula
%% New version
% Twist the XY plane about z = [0, 0, 1] by phi = atan2(vy, vx)
% using tangent vector v = [vx, vy], then tilt from z axis to normal
% vector n = [nx, ny, nz] so that the circle's normal is n regardless
% of v = [vx, vy]
v = [T.vx, T.vy];
n = [T.nx, T.ny, T.nz];
R = calcrotation(T,v,n);

% Translate to center c = [cx, cy, cz]
c = [T.cx, T.cy, T.cz];

for i=1:8
    P{i} = P{i}*R + c;
end
