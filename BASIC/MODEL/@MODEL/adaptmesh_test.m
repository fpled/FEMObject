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

G = appendnodedata(G,q);

options = '';
if ischarin('gmshoptions',varargin)
    options = [options ' ' getcharin('gmshoptions',varargin)];
end

if isempty(strfind(options,'-bmg')) % for compatibility with Matlab version < 9.1 (R2016b)
% if ~contains(options,'-bmg') % for Matlab versions >= 9.1 (R2016b)
    options = [options ' -bmg ' getfile(G,'.msh') '.pos'];
end

G = remesh(G,dim,'gmshoptions',options);

n = max(nargout,1);
varargout = cell(1,n);
[varargout{:}] = gmsh2femobject(indim,G,dim:-1:dim-n+1,varargin{:});
