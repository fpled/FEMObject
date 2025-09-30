function u = compoundspline(u,numbercurve,values,numIntervals)
% function u = compoundspline(u,numbercurve,values,numIntervals)

if nargin<4 || isempty(numIntervals)
    numIntervals = 5;
end

u = compoundusing(u,'Spline',numbercurve,values,numIntervals);
