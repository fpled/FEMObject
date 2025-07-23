function varargout = gmshLshape(cl,filename,indim,varargin)
% function varargout = gmshLshape(cl,filename,indim)
% cl : characteristic length
% filename : file name (optional)
% indim : space dimension (optional, 2 by default)

if nargin<3 || isempty(indim)
    indim = 2;
end

P{1} = [0, 0, 0];
P{2} = [1, 0, 0];
P{3} = [1, 1, 0];
P{4} = [2, 1, 0];
P{5} = [2, 2, 0];
P{6} = [0, 2, 0];

G = GMSHFILE();
if nargin>=2 && ischar(filename)
    G = setfile(G,filename);
end

G = createpoints(G,P,cl,1:6);
G = createcontour(G,1:6,1:6,7);
G = createplanesurface(G,7,1);
if indim==3
    vect = [0,0,1];
    [G,out] = extrude(G,vect,'Surface',1,varargin{:});
    numbervolume = [out,'[1]'];
end
if ischarin('recombine',varargin)
    if indim==2
        G = recombinesurface(G,1);
    elseif indim==3
        G = recombinesurface(G);
    end
end
if indim==2
    G = createphysicalsurface(G,1,1);
elseif indim==3
    G = createphysicalvolume(G,numbervolume,1);
end

varargin = delonlycharin('recombine',varargin);

n = max(nargout,1);
varargout = cell(1,n);
[varargout{:}] = gmsh2femobject(indim,G,indim:-1:indim-n+1,varargin{:});
