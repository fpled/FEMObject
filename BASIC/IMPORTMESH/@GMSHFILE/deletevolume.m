function u = deletevolume(u,numbervolume,varargin)
% function u = deletevolume(u,numbervolume,varargin)

u = delete(u,'Volume',numbervolume,varargin{:});
