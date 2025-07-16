function u = createellipsecontour(u,numbercenter,numberpoints,numberlines,numberlineloop,varargin)
% function u = createellipsecontour(u,numcenter,numberpoints,numberlines,numberlineloop,reverse)

if nargin<6
    reverse = 1;
else
    reverse = getcharin('reverse',varargin,1);
end

seg = [1:length(numberpoints);2:length(numberpoints),1];
maj = numberpoints;
seg = numberpoints(seg)';

u = createellipses(u,numbercenter,seg,maj,numberlines);
u = createcurveloop(u,reverse*numberlines,numberlineloop);
