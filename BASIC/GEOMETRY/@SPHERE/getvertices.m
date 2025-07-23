function P = getvertices(S)
% function P = getvertices(S)

% Radius
r = S.r;

P = cell(1,6);
% Vertices
P{1} = [ r,  0,  0]; % +x
P{2} = [ 0,  r,  0]; % +y
P{3} = [ 0,  0,  r]; % +z
P{4} = [-r,  0,  0]; % -x
P{5} = [ 0, -r,  0]; % -y
P{6} = [ 0,  0, -r]; % -z

%% Old version
% Rotate around axis n = [nx, ny, nz] by angle of rotation phi =
% atan2(vy, vx) using tangent vector v = [vx, vy] via Rodrigues'
% rotation formula
%% New version
% Twist the XY plane about z = [0, 0, 1] by phi = atan2(vy, vx)
% using tangent vector v = [vx, vy], then tilt from z axis to normal
% vector n = [nx, ny, nz] so that the circle's normal is n regardless
% of v = [vx, vy]
v = [S.vx, S.vy];
n = [S.nx, S.ny, S.nz];
R = calcrotation(S,v,n);

% Translate to center c = [cx, cy, cz]
c = [S.cx, S.cy, S.cz];

for i=1:6
    P{i} = P{i}*R + c;
end