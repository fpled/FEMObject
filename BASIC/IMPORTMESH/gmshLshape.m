function varargout = gmshLshape(cl,filename,indim,varargin)
% function varargout = gmshLshape(cl,filename,indim)
% cl : characteristic length
% filename : file name (optional)
% indim : space dimension (optional, 2 by default)

if nargin<3 || isempty(indim)
    indim = 2;
end

dim = indim;

P{1} = [0, 0, 0];
P{2} = [1, 0, 0];
P{3} = [1, 1, 0];
P{4} = [2, 1, 0];
P{5} = [2, 2, 0];
P{6} = [0, 2, 0];

br = @(tag,k) sprintf('%s[%d]', tag, k); % bracket reference helper

G = GMSHFILE();
if nargin>=2 && ischar(filename)
    G = setfile(G,filename);
end

numpoints = 1:6;
numlines = 1:6;
numlineloop = 1;
numsurface = 1;
G = createpoints(G,P,cl,numpoints);
G = createcontour(G,numpoints,numlines,numlineloop);
G = createplanesurface(G,numlineloop,numsurface);
if indim==3
    t = 1;
    vect = [0,0,t];
    if ischarin('recombine',varargin) && ~ischarin('Layers',varargin)
        numlayers = max(1, round(t/cl));
        varargin = [varargin, {'Layers',numlayers}];
    end
    [G,tag] = extrude(G,vect,'Surface',numsurface,varargin{:});
    numtopsurface = br(tag,0); % top surface
    numvolume     = br(tag,1); % volume
end
if ischarin('recombine',varargin)
    if indim==2
        G = recombinesurface(G,numsurface);
    elseif indim==3
        G = recombinesurface(G,{numsurface,numtopsurface});
    end
end
if indim==2
    numphysicalsurface = 1;
    G = createphysicalsurface(G,numsurface,numphysicalsurface);
elseif indim==3
    numphysicalvolume = 1;
    G = createphysicalvolume(G,numvolume,numphysicalvolume);
end

varargin = delonlycharin('recombine',varargin);

n = max(nargout,1);
varargout = cell(1,n);
[varargout{:}] = gmsh2femobject(indim,G,dim:-1:dim-n+1,varargin{:});
