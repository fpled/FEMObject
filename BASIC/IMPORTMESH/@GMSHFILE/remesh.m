function u = remesh(u,options,varargin)
% function u = remesh(u,dim)
% launch gmsh u.msh -dim

% function u = remesh(u,dim,'gmshoptions',options)
% options: char containing gmsh options
% launch gmsh u.msh -dim options

% function u = remesh(u,[],'gmshoptions',options)
% options: char containing gmsh options
% launch gmsh u.msh options

if nargin==1
    error('rentrer les options de gsmh')
end

if verLessThan('matlab','9.1') % compatibility (<R2016b)
    contain = @(str,pat) ~isempty(strfind(str,pat));
else
    contain = @contains;
end

if ~isempty(options) && isa(options,'double')
    options = [' -' num2str(options)];
else
    options = '';
end

if ischarin('gmshoptions',varargin)
    options = [options ' ' getcharin('gmshoptions',varargin)];
end

if ~contain(options,'-format')
    options = [options ' -format msh2'];
end

filegeo = getfile(u,'.geo','post');
if exist(filegeo,'file')
    options = [filegeo options];
end

u = runfile(u,'.msh',options);

fileopt = getfile(u,'.msh.opt');
if exist(fileopt,'file')
    u = deletefile(u,'.msh.opt');
end
