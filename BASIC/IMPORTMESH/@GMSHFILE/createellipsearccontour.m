function u = createellipsearccontour(u,numbercenter,numberpoints,numbermajorpoints,numbercurves,numbercurveloop,varargin)
% function u = createellipsearccontour(u,numcenter,numberpoints,numbermajorpoints,numbercurves,numbercurveloop,reverse)

if nargin<7
    reverse = 1;
else
    reverse = getcharin('reverse',varargin,1);
end

% Open arc
n = length(numberpoints);
seg_arc = [1:n-1; 2:n]';
seg_arc = numberpoints(seg_arc);

% Radial lines
seg_line = [numbercenter, numberpoints(1); ...
            numberpoints(end), numbercenter];

numberlines = numbercurves([1,end]);
numberarcs  = numbercurves(2:end-1);

if length(numbermajorpoints)~=length(numberarcs)
    error('numbermajorpoints must be the same length as numberarcs.');
end

u = createlines(u,seg_line,numberlines); % radial lines
u = createellipsearcs(u,numbercenter,numbermajorpoints,seg_arc,numberarcs); % middle arcs
if reverse==-1
    numbercurves = fliplr(numbercurves);
end
u = createcurveloop(u,numbercurves,numbercurveloop);
