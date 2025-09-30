function u = deletepoint(u,numberpoint,varargin)
% function u = deletepoint(u,numberpoint,varargin)

u = delete(u,'Point',numberpoint,varargin{:});
