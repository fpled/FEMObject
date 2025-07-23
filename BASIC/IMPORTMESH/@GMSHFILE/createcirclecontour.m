function u = createcirclecontour(u,numbercenter,numberpoints,numbercurves,numbercurveloop,varargin)
% function u = createcirclecontour(u,numcenter,numberpoints,numbercurves,numbercurveloop,reverse)

if nargin<6
    reverse = 1;
else
    reverse = getcharin('reverse',varargin,1);
end

% Closed loop
n = length(numberpoints);
seg = [1:n; 2:n,1]';
seg = numberpoints(seg);

u = createcirclearcs(u,numbercenter,seg,numbercurves);
u = createcurveloop(u,reverse*numbercurves,numbercurveloop);
