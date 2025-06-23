function u = remesh(u,options,varargin)
% function u = remesh(u,dim)
% launch gmsh u -dim

% function u = remesh(u,[],'gmshoptions',options)
% options: char
% launch gmsh u options

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

u = runfile(u,'.msh',options);
