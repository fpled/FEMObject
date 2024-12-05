function varargout = gmshdomainwithinteriorcrack(D,C,clD,clC,filename,indim,varargin)
% function varargout = gmshdomainwithinteriorcrack(D,C,clD,clC,filename,indim)
% D : DOMAIN
% C : LIGNE in dim 2, QUADRANGLE in dim 3
% clD, clC : characteristic lengths
% filename : file name (optional)
% indim : space dimension (optional, getindim(D) by default)

noduplicate = ischarin('noduplicate',varargin);
varargin = delonlycharin('noduplicate',varargin);

if nargin<6 || isempty(indim)
    indim = getindim(D);
end
if nargin<4 || isempty(clC)
    clC = clD;
end

if ~iscell(C)
    C = {C};
end
if length(clC)==1
    clC = repmat(clC,1,length(C));
end

if indim==2
    G = gmshfile(D,clD,1:4,1:4,1,1,varargin{:});
    numpoints = 4+(1:2*length(C));
    numlines = 4+(1:length(C));
    for j=1:length(C)
        GC = gmshfile(C{j},clC(j),numpoints([2*j-1,2*j]),numlines(j));
        G = G+GC;
        G = embedlineinsurface(G,numlines(j),1);
    end
    if ~noduplicate
        physicalgroup = 1;
        G = createphysicalpoint(G,numpoints,1);
        G = createphysicalline(G,numlines,physicalgroup);
        G = createphysicalsurface(G,1,1);
    end
elseif indim==3
    G = gmshfile(D,clD,1:8,1:4,1,1,varargin{:});
    numpoints = 8+(1:4*length(C));
    numlines = 4+(1:4*length(C));
    numbersurface = 1+(1:length(C));
    for j=1:length(C)    
        GC = gmshfile(C{j},clC(j),numpoints(4*j-3:4*j),numlines(4*j-3:4*j),j+1,numbersurface(j),varargin{:});
        G = G+GC;
        G = embedsurfaceinvolume(G,numbersurface(j),1);
    end
    if ~noduplicate
        physicalgroup = 1;
        G = createphysicalline(G,numlines,1);
        G = createphysicalsurface(G,numbersurface,physicalgroup);
        G = createphysicalvolume(G,1,1);
    end
end
varargin = delonlycharin('recombine',varargin);

if nargin>=5 && ischar(filename)
    G = setfile(G,filename);
end

n=max(nargout,1);
varargout = cell(1,n);
[varargout{:}] = gmsh2femobject(indim,G,getdim(D):-1:getdim(D)-n+1,varargin{:});

if ~noduplicate
    G = createcrack(G,getdim(D)-1,physicalgroup);
    G = remesh(G,getdim(D),varargin{:});
    G = deleteoptfile(G);
    
    [varargout{:}] = gmsh2femobject(indim,G,getdim(D):-1:getdim(D)-n+1,varargin{:});
end