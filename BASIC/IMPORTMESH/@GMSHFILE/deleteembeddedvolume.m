function u = deleteembeddedvolume(u,numbervolume,varargin)
% function u = deleteembeddedvolume(u,numbervolume,varargin)

u = deleteembedded(u,'Volume',numbervolume,varargin{:});
