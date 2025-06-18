function u = createellipsecontour(u,numbercenter,numberpoints,numberlines,numberlineloop)
% function u = createellipsecontour(u,numcenter,numberpoints,numberlines,numberlineloop)

seg = [1:length(numberpoints);2:length(numberpoints),1];
maj = numberpoints;
seg = numberpoints(seg)';

u = createellipses(u,numbercenter,seg,maj,numberlines);
u = createcurveloop(u,numberlines,numberlineloop);
