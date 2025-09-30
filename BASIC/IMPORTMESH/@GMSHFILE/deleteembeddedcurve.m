function u = deleteembeddedcurve(u,numbercurve,varargin)
% function u = deleteembeddedcurve(u,numbercurve,varargin)

u = deleteembedded(u,'Curve',numbercurve,varargin{:});
