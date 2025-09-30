function u = compoundBspline(u,numbercurve,values,numIntervals)
% function u = compoundBspline(u,numbercurve,values,numIntervals)

if nargin<4 || isempty(numIntervals)
    numIntervals = 20;
end

u = compoundusing(u,'BSpline',numbercurve,values,numIntervals);
