function u = deleteembeddedsurface(u,numbersurface,varargin)
% function u = deleteembeddedsurface(u,numbersurface,varargin)

u = deleteembedded(u,'Surface',numbersurface,varargin{:});
