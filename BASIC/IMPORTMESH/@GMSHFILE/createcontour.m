function u = createcontour(u,numberpoints,numberlines,numberlineloop)
% function u = createcontour(u,numberpoints,numberlines,numberlineloop)

seg = [1:length(numberpoints);2:length(numberpoints),1]';
seg = numberpoints(seg);
u = createlines(u,seg,numberlines);
u = createcurveloop(u,numberlines,numberlineloop);
