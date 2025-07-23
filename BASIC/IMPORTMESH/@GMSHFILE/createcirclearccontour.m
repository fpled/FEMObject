function u = createcirclearccontour(u,numbercenter,numberpoints,numbercurves,numbercurveloop,varargin)
% function u = createcirclearccontour(u,numcenter,numberpoints,numbercurves,numbercurveloop,reverse)

if nargin<6
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

u = createlines(u,seg_line,numberlines); % radial lines
u = createcirclearcs(u,numbercenter,seg_arc,numberarcs); % middle arcs
if reverse==-1
    numbercurves = fliplr(numbercurves);
end
u = createcurveloop(u,numbercurves,numbercurveloop);
