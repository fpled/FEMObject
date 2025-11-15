function u = mesh(u,options,varargin)
% function u = mesh(u,dim)
% launch gmsh u.geo -dim

% function u = mesh(u,dim,'gmshoptions',options)
% options: char containing gmsh options
% launch gmsh u.geo -dim options

% function u = mesh(u,[],'gmshoptions',options)
% options: char containing gmsh options
% launch gmsh u.geo options

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

u = runfile(u,'.geo',options);
