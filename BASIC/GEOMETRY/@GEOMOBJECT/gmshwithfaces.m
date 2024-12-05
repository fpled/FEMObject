function varargout = gmshwithfaces(D,varargin)
% function G = gmshwithfaces(D,varargin)
% D : GEOMOBJECT

filename = getcharin('filename',varargin,'gmsh_file');
indim = getcharin('indim',varargin,getindim(D));
recombine = ischarin('recombine',varargin);
varargin = delcharin({'filename','indim'},varargin);
varargin = delonlycharin('recombine',varargin);

[G,numbersurface] = gmshfile(D,varargin{:});
if recombine
    G = recombinesurface(G,numbersurface);
end
G = setfile(G,filename);
n=max(nargout,1);
varargout=cell(1,n);
[varargout{:}] = gmsh2femobject(indim,G,getdim(D):-1:getdim(D)-1);