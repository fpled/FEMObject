function u = createcirclecontour(u,numbercenter,numberpoints,numberlines,numberlineloop,varargin)
% function u = createcirclecontour(u,numcenter,numberpoints,numberlines,numberlineloop,reverse)

if nargin<6
    reverse=1;
else
    reverse = getcharin('reverse',varargin,1);
end

seg = [1:length(numberpoints);2:length(numberpoints),1];
seg = numberpoints(seg)';

u = createcircles(u,numbercenter,seg,numberlines);
u = createcurveloop(u,reverse*numberlines,numberlineloop);
