function varargout = adaptmesh(M,q,filename,varargin)
% function varargout = adaptmesh(M,q,filename,varargin)

if ~isa(M,'MODEL')
    varargout = cell(1,nargout);
    [varargout{:}] = adaptmesh(q,M,filename,varargin{:});
    return
end

if israndom(q)
    error('la solution est aleatoire')
end
q = unfreevector(M,q);

indim = getindim(M);
dim = getdim(M);

G = GMSHFILE();
if nargin>=3 && ischar(filename)
    G = setfile(G,filename);
end

if ischarin('mmgoptions',varargin)
    % Mesh adaptation using Mmg
    G = adaptmesh_mmg(G,q,dim,varargin{:});
else
    % Mesh adaptation using Gmsh
    G = adaptmesh_gmsh(G,q,dim,varargin{:});
end

n = max(nargout,1);
varargout = cell(1,n);
[varargout{:}] = gmsh2femobject(indim,G,dim:-1:dim-n+1,varargin{:});
