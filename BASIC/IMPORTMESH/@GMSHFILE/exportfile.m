function u = exportfile(u,ext,format,varargin)
% function u = exportfile(u)
% launch gmsh u.msh -save -format msh2

% function u = exportfile(u,ext)
% ext: char containing file extension
% launch gmsh [u ext] -save -format msh2

% function u = exportfile(u,ext,format)
% format: char containing file format
% launch gmsh [u ext] -save -format format

% function u = exportfile(u,ext,[],'gmshoptions',options)
% options: char containing gmsh options
% launch gmsh [u ext] -save options

% function u = exportfilemsh(u,ext,format,'gmshoptions',options)
% options: char containing gmsh options
% launch gmsh [u ext] -save -format format options

if nargin<2 || isempty(ext)
    ext = '.msh';
end
if nargin<3 || isempty(format)
    format = 'msh2';
end

options = [' -save -format ' format];
% options = [' -save -format ' format ' -preserve_numbering_msh2'];

if ischarin('gmshoptions',varargin)
    options = [options ' ' getcharin('gmshoptions',varargin)];
end

u = runfile(u,ext,options);
