function P = getvertices(E)
% function P = getvertices(E)

% Semi-axes
a = E.a;
b = E.b;
c = E.c;

% Vertices
P = cell(1,6);
P{1} = [ a,  0,  0]; % +x
P{2} = [ 0,  b,  0]; % +y
P{3} = [ 0,  0,  c]; % +z
P{4} = [-a,  0,  0]; % -x
P{5} = [ 0, -b,  0]; % -y
P{6} = [ 0,  0, -c]; % -z

% Rotation matrix
v = [E.vx, E.vy];
n = [E.nx, E.ny, E.nz];
R = calcrotation(E,v,n);

% Center
c = [E.cx, E.cy, E.cz];

% Rotate into global frame and translate to center
for i=1:6
    P{i} = P{i}*R + c;
end
