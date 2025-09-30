function [u,tag] = symmetrize(u,coefficients,name,number,varargin)
% function u = symmetrize(u,coefficients,name,number[,tag,'duplicate'])
% function [u,tag] = symmetrize(u,coefficients,name,number[,tag,'duplicate'])
% tag = 'out' by default if nargout==2
% This is a deprecated synonym for mirror.

[u,tag] = mirror(u,coefficients,name,number,varargin{:});
