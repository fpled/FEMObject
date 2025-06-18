function varargout = exportmesh(M,filename,ext,format,varargin)
% function varargout = exportmesh(M,filename,ext,format,varargin)

indim = getindim(M);
dim = getdim(M);

G = GMSHFILE();
if nargin>=2 && ischar(filename)
    G = setfile(G,filename);
end

if nargin<3
    ext = '.msh';
end
if nargin<4
    format = 'msh2';
end

G = exportfile(G,ext,format,varargin{:});

n=max(nargout,1);
varargout = cell(1,n);
[varargout{:}] = gmsh2femobject(indim,G,dim:-1:dim-n+1,varargin{:});
