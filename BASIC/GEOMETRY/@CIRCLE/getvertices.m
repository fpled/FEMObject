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
        
        % Rotate in-plane around z = [0, 0, 1] by angle of rotation theta =
        % atan2(vy, vx) using tangent vector v = [vx, vy]
        v = [C.vx, C.vy];
        R = calcrotation(C,v);
        
        % Translate to center c = [cx, cy]
        c = [C.cx, C.cy];
        
    case 3
        % Vertices
        P{1} = [ r,  0, 0]; % +x
        P{2} = [ 0,  r, 0]; % +y
        P{3} = [-r,  0, 0]; % -x
        P{4} = [ 0, -r, 0]; % -y
        
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
        
    otherwise
        error('Wrong space dimension');
end

for i=1:4
    P{i} = P{i}*R + c;
end
