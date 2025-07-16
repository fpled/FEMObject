function P = getvertices(E)
% function P = getvertices(E)

% Semi-axes
a = E.a;
b = E.b;
c = E.c;

P = cell(1,6);
% Vertices
P{1} = [-a,  0,  0];
P{2} = [ 0, -b,  0];
P{3} = [ 0,  0, -c];
P{4} = [ a,  0,  0];
P{5} = [ 0,  b,  0];
P{6} = [ 0,  b,  c];

%% Old version
% Rotate around axis n = [nx, ny, nz] by angle of rotation phi =
% atan2(vy, vx) using tangent vector v = [vx, vy] via Rodrigues'
% rotation formula
%% New version
% Twist the XY plane about z = [0, 0, 1] by phi = atan2(vy, vx)
% using tangent vector v = [vx, vy], then tilt from z axis to normal
% vector n = [nx, ny, nz] so that the circle's normal is n regardless
% of v = [vx, vy]
v = [E.vx, E.vy];
n = [E.nx, E.ny, E.nz];
R = calcrotation(E,v,n);

% Translate to center c = [cx, cy, cz]
c = [E.cx, E.cy, E.cz];

for i=1:6
    P{i} = P{i}*R + c;
end
