function P = getvertices(S)
% function P = getvertices(S)

% Radius
r = S.r;

% Vertices
P = cell(1,6);
P{1} = [ r,  0,  0]; % +x
P{2} = [ 0,  r,  0]; % +y
P{3} = [ 0,  0,  r]; % +z
P{4} = [-r,  0,  0]; % -x
P{5} = [ 0, -r,  0]; % -y
P{6} = [ 0,  0, -r]; % -z

% Rotation matrix
v = [S.vx, S.vy];
n = [S.nx, S.ny, S.nz];
R = calcrotation(S,v,n);

% Center
c = [S.cx, S.cy, S.cz];

% Rotate into global frame and translate to center
for i=1:6
    P{i} = P{i}*R + c;
end