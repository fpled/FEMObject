function [rep,P] = ispointin(C,P)
% function [rep,P] = ispointin(C,P)

c = getcoord(P);
 
% Radius
r = C.r;

% Spherical condition: point lies within sphere radius
d = sqrt((c(:,1) - C.cx).^2 + (c(:,2) - C.cy).^2 + (c(:,3) - C.cz).^2); % Euclidean distance to center
inSphere = d <= r + eps;
rep = find(inSphere);

if nargout==2
    P = P(rep);
end
