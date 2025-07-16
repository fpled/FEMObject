function [rep,P] = ispointin(E,P)
% function [rep,P] = ispointin(E,P)

c = getcoord(P);

% Semi-axes
a = E.a;
b = E.b;

switch E.indim
    case 2
        % Radial condition: point lies within ellipse semi-axes
        xc = (c(:,1)-E.cx)/a;
        yc = (c(:,2)-E.cy)/b;
        inEllipse = xc.^2 + yc.^2 <= 1 + eps;
        rep = find(inEllipse);
        
    case 3
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
        
        % Apply inverse transform to point: project into the ellipse local frame
        center = [E.cx, E.cy, E.cz];
        vec = c - center; % vector from center to point
        vec = vec * R';
        
        % In-plane condition: point lies in the local z=0 plane
        tol = getfemobjectoptions('tolerancepoint');
        inPlane = abs(vec(:,3)) < tol;
        
        % Radial condition: point lies within ellipse semi-axes
        xc = vec(:,1)/a;
        yc = vec(:,2)/b;
        inEllipse = (xc.^2 + yc.^2) <= 1 + eps;
        
        % Both conditions: point lies in the plane and within ellipse semi-axes
        rep = find(inPlane & inEllipse);
        
    otherwise
        error('Wrong space dimension');
end

if nargout==2
    P = P(rep);
end
