function u = createcontour(u,numberpoints,numberlines,numberlineloop)
% function u = createcontour(u,numberpoints,numberlines,numberlineloop)

% Closed loop
n = length(numberpoints);
seg = [1:n; 2:n,1]';
seg = numberpoints(seg);

u = createlines(u,seg,numberlines);
u = createcurveloop(u,numberlines,numberlineloop);
