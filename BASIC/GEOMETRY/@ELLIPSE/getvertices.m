function P = getvertices(E)
% function P = getvertices(E)

% Semi-axes
a = E.a;
b = E.b;

P = cell(1,4);
switch E.indim
    case 2
        % Vertices
        P{1} = [ a,  0]; % +x
        P{2} = [ 0,  b]; % +y
        P{3} = [-a,  0]; % -x
        P{4} = [ 0, -b]; % -y
        
        % Rotation matrix
        v = [E.vx, E.vy];
        R = calcrotation(E,v);
        
        % Center
        c = [E.cx, E.cy];
        
    case 3
        % Vertices
        P{1} = [ a,  0, 0]; % +x
        P{2} = [ 0,  b, 0]; % +y
        P{3} = [-a,  0, 0]; % -x
        P{4} = [ 0, -b, 0]; % -y
        
        % Rotation matrix
        v = [E.vx, E.vy];
        n = [E.nx, E.ny, E.nz];
        R = calcrotation(E,v,n);
        
        % Center
        c = [E.cx, E.cy, E.cz];
        
    otherwise
        error('Wrong space dimension');
end

% Rotate into global frame and translate to center
for i=1:4
    P{i} = P{i}*R + c;
end
