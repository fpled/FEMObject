function u = deleteline(u,numberline,varargin)
% function u = deleteline(u,numberline,varargin)

u = delete(u,'Line',numberline,varargin{:});
