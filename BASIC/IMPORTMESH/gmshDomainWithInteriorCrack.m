function varargout = gmshDomainWithInteriorCrack(D,C,clD,clC,filename,indim,varargin)
% function varargout = gmshDomainWithInteriorCrack(D,C,clD,clC,filename,indim)
% D : DOMAIN
% C : LINE in dim 2, QUADRANGLE in dim 3
% clD, clC : characteristic lengths
% filename : file name (optional)
% indim : space dimension (optiownal, getindim(D) by default)

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
if isscalar(clC)
    clC = repmat(clC,1,length(C));
end

if indim==2
    numpoints = 1:4;
    numlines = 1:4;
    numlineloop = 1;
    numsurface = 1;
    G = gmshfile(D,clD,numpoints,numlines,numlineloop,numsurface,varargin{:});
    % numembeddedpoints = [];
    numembeddedcurves = [];
    for j=1:length(C)
        numpoints = numpoints(end)+(1:2);
        numlines = numlines(end)+1;
        GC = gmshfile(C{j},clC(j),numpoints,numlines);
        G = G+GC;
        G = embedcurveinsurface(G,numlines,numsurface);
        % numembeddedpoints = [numembeddedpoints,numpoints];
        numembeddedcurves = [numembeddedcurves,numlines];
    end
    
    if ~noduplicate
        physicalgroup = 1;
        % G = createphysicalpoint(G,numembeddedpoints,1);
        G = createphysicalcurve(G,numembeddedcurves,physicalgroup);
    end
    
    numphysicalsurface = 1;
    G = createphysicalsurface(G,numsurface,numphysicalsurface);

elseif indim==3
    numpoints = 1:8;
    numlines = 1:12;
    numlineloop = 1:6;
    numsurface = 1:6;
    numsurfaceloop = 1;
    numvolume = 1;
    G = gmshfile(D,clD,numpoints,numlines,numlineloop,numsurface,numsurfaceloop,numvolume,varargin{:});
    % numembeddedpoints = [];
    % numembeddedcurves = [];
    numembeddedsurfaces = [];
    for j=1:length(C)
        numpoints = numpoints(end)+(1:4);
        numlines = numlines(end)+(1:4);
        numlineloop = numsurfaceloop(end)+1;
        numsurface = numsurface(end)+1;
        GC = gmshfile(C{j},clC(j),numpoints,numlines,numlineloop,numsurface,varargin{:});
        G = G+GC;
        G = embedsurfaceinvolume(G,numsurface,numvolume);
        % numembeddedpoints = [numembeddedpoints,numpoints];
        % numembeddedcurves = [numembeddedcurves,numlines];
        numembeddedsurfaces = [numembeddedsurfaces,numsurface];
    end
    
    if ~noduplicate
        physicalgroup = 1;
        % G = createphysicalpoint(G,numembeddedpoints,1);
        % G = createphysicalcurve(G,numembeddedcurves,1);
        G = createphysicalsurface(G,numembeddedsurfaces,physicalgroup);
    end
    
    numphysicalvolume = 1;
    G = createphysicalvolume(G,numvolume,numphysicalvolume);

end

varargin = delonlycharin('recombine',varargin);

if nargin>=5 && ischar(filename)
    G = setfile(G,filename);
end

n = max(nargout,1);
varargout = cell(1,n);
[varargout{:}] = gmsh2femobject(indim,G,getdim(D):-1:getdim(D)-n+1,varargin{:});

if ~noduplicate
    G = createcrack(G,getdim(D)-1,physicalgroup);
    G = remesh(G,getdim(D),varargin{:});
    G = deleteoptfile(G);
    
    [varargout{:}] = gmsh2femobject(indim,G,getdim(D):-1:getdim(D)-n+1,varargin{:});
end
