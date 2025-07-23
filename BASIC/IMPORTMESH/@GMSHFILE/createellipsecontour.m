function u = createellipsecontour(u,numbercenter,numberpoints,numbermajorpoints,numbercurves,numbercurveloop,varargin)
% function u = createellipsecontour(u,numcenter,numberpoints,numbermajorpoints,numbercurves,numbercurveloop,reverse)

if nargin<7
    reverse = 1;
else
    reverse = getcharin('reverse',varargin,1);
end

if length(numbermajorpoints)~=length(numbercurves)
    error('numbermajorpoints must be the same length as numbercurves.');
end

% Closed loop
n = length(numberpoints);
seg = [1:n; 2:n,1]';
seg = numberpoints(seg);

u = createellipsearcs(u,numbercenter,seg,numbermajorpoints,numbercurves);
u = createcurveloop(u,reverse*numbercurves,numbercurveloop);
