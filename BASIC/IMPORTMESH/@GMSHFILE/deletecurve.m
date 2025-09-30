function u = deletecurve(u,numbercurve,varargin)
% function u = deletecurve(u,numbercurve,varargin)

u = delete(u,'Curve',numbercurve,varargin{:});
