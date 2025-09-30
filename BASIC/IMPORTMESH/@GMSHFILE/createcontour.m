function u = createcontour(u,numberpoints,numberlines,numberlineloop,varargin)
% function u = createcontour(u,numberpoints,numberlines,numberlineloop)

if nargin<4, numberlineloop = []; end

reverse = getcharin('reverse',varargin,1);
assert(ismember(reverse,[-1,1]),'''reverse'' must be +1 or -1');

% Closed loop
n = length(numberpoints);
seg = [1:n; 2:n,1]';
seg = numberpoints(seg);

u = createlines(u,seg,numberlines);

if ~isempty(numberlineloop)
    numberlines = reverse * numberlines;
    u = createcurveloop(u,numberlines,numberlineloop);
end
