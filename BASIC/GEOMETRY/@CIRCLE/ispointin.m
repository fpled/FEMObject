function [rep,P] = ispointin(C,P)
% function [rep,P] = ispointin(C,P)

c = getcoord(P);

% Radius
r = C.r;

switch C.indim
    case 2
        % Radial condition: point lies within circle radius
        % d = sqrt((c(:,1)-C.cx).^2 + (c(:,2)-C.cy).^2); % Euclidian distance to center
        d = hypot(c(:,1)-C.cx, c(:,2)-C.cy);           % Euclidian distance to center
        inCircle = d <= r + eps;
        rep = find(inCircle);
        
    case 3
        % Rotation matrix
        v = [C.vx, C.vy];
        n = [C.nx, C.ny, C.nz];
        R = calcrotation(C,v,n);
        
        % Apply inverse transform to point: project into the circle local frame
        center = [C.cx, C.cy, C.cz];
        vec = c - center; % vector from center to point
        vec = vec * R';
        
        % In-plane condition: point lies in the local z=0 plane
        tol = getfemobjectoptions('tolerancepoint');
        inPlane = abs(vec(:,3)) < tol;
        
        % Radial condition: point lies within circle radius
        % d = sqrt(vec(:,1).^2 + vec(:,2).^2); % radial distance in the plane
        d = hypot(vec(:,1), vec(:,2));       % radial distance in the plane
        inCircle = d <= r + eps;
        
        % Both conditions: point lies in the plane and within circle radius
        rep = find(inPlane & inCircle);
        
    otherwise
        error('Wrong space dimension');
end

if nargout==2
    P = P(rep);
end
