function u = createellipsecontour(u,numbercenter,numberpoints,numbermajorpoints,numbercurves,numbercurveloop,varargin)
% function u = createellipsecontour(u,numcenter,numberpoints,numbermajorpoints,numbercurves,numbercurveloop,reverse)

if nargin<6, numbercurveloop = []; end

reverse = getcharin('reverse',varargin,1);
assert(ismember(reverse,[-1,1]),'''reverse'' must be +1 or -1');

if length(numbermajorpoints)~=length(numbercurves)
    error('numbermajorpoints must be the same length as numbercurves.');
end

% Closed loop
n = length(numberpoints);
seg = [1:n; 2:n,1]';
seg = numberpoints(seg);

u = createellipsearcs(u,numbercenter,seg,numbermajorpoints,numbercurves);

if ~isempty(numbercurveloop)
    numbercurves = reverse * numbercurves;
    u = createcurveloop(u,numbercurves,numbercurveloop);
end
