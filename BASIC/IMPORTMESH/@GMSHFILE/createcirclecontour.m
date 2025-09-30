function u = createcirclecontour(u,numbercenter,numberpoints,numbercurves,numbercurveloop,varargin)
% function u = createcirclecontour(u,numcenter,numberpoints,numbercurves,numbercurveloop,reverse)

if nargin<5, numbercurveloop = []; end

reverse = getcharin('reverse',varargin,1);
assert(ismember(reverse,[-1,1]),'''reverse'' must be +1 or -1');

% Closed loop
n = length(numberpoints);
seg = [1:n; 2:n,1]';
seg = numberpoints(seg);

u = createcirclearcs(u,numbercenter,seg,numbercurves);

if ~isempty(numbercurveloop)
    numbercurves = reverse * numbercurves;
    u = createcurveloop(u,numbercurves,numbercurveloop);
end
