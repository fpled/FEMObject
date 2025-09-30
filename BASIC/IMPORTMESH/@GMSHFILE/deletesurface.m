function u = deletesurface(u,numbersurface,varargin)
% function u = deletesurface(u,numbersurface,varargin)

u = delete(u,'Surface',numbersurface,varargin{:});
