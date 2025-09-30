function P = getvertices(C)
% function P = getvertices(C)

% Radius
r = C.r;

P = cell(1,4);
switch C.indim
    case 2
        % Vertices
        P{1} = [ r,  0]; % +x
        P{2} = [ 0,  r]; % +y
        P{3} = [-r,  0]; % -x
        P{4} = [ 0, -r]; % -y
        
        % Rotation matrix
        v = [C.vx, C.vy];
        R = calcrotation(C,v);
        
        % Center
        c = [C.cx, C.cy];
        
    case 3
        % Vertices
        P{1} = [ r,  0, 0]; % +x
        P{2} = [ 0,  r, 0]; % +y
        P{3} = [-r,  0, 0]; % -x
        P{4} = [ 0, -r, 0]; % -y
        
        % Rotation matrix
        v = [C.vx, C.vy];
        n = [C.nx, C.ny, C.nz];
        R = calcrotation(C,v,n);
        
        % Center
        c = [C.cx, C.cy, C.cz];
        
    otherwise
        error('Wrong space dimension');
end

% Rotate into global frame and translate to center
for i=1:4
    P{i} = P{i}*R + c;
end
