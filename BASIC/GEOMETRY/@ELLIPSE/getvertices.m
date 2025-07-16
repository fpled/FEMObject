function P = getvertices(E)
% function P = getvertices(E)

% Semi-axes
a = E.a;
b = E.b;

P = cell(1,4);
switch E.indim
    case 2
        % Vertices
        P{1} = [-a,  0];
        P{2} = [ 0, -b];
        P{3} = [ a,  0];
        P{4} = [ 0,  b];
        
        % Rotate in-plane around z = [0, 0, 1] by angle of rotation theta =
        % atan2(vy, vx) using tangent vector v = [vx, vy]
        v = [E.vx, E.vy];
        R = calcrotation(E,v);
        
        % Translate to center c = [cx, cy]
        c = [E.cx, E.cy];
        
    case 3
        % Vertices
        P{1} = [-a,  0, 0];
        P{2} = [ 0, -b, 0];
        P{3} = [ a,  0, 0];
        P{4} = [ 0,  b, 0];
        
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
        
    otherwise
        error('Wrong space dimension');
end

for i=1:4
    P{i} = P{i}*R + c;
end
