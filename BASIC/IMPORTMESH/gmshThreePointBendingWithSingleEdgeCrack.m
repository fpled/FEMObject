function varargout = gmshThreePointBendingWithSingleEdgeCrack(L,H,ls,w,a,b,t,clD,clC,clS,filename,indim,varargin)
% function varargout = gmshThreePointBendingWithSingleEdgeCrack(L,H,ls,w,a,b,t,clD,clC,clS,filename,indim)
% L : length
% H : height
% ls : location of the supports from the outer edges
% w : flat punch width
% a : length of the edge crack (or notch)
% b : location of the edge crack (or notch) from the left outer edge
% t : thickness
% clD, clC, clS : characteristic lengths
% filename : file name (optional)
% indim : space dimension (optional, 2 by default)

Box         = getcharin('Box',varargin,[]);
noduplicate = ischarin('noduplicate',varargin);
refinecrack = ischarin('refinecrack',varargin);
recombine   = ischarin('recombine',varargin);

varargin = delcharin('Box',varargin);
varargin = delonlycharin({'noduplicate','refinecrack','recombine'},varargin);

if nargin<12 || isempty(indim)
    indim = 2;
end
if nargin<10 || isempty(clC)
    clS = clD;
end
if nargin<9 || isempty(clC)
    clC = clD;
end
if nargin<7 || isempty(t)
    t = 1;
end

dim = 2;

% P{1}  = [0       , 0];
% P{2}  = [ls-w/2  , 0];
% P{3}  = [ls      , 0];
% P{4}  = [ls+w/2  , 0];
% P{5}  = [L-ls-w/2, 0];
% P{6}  = [L-ls    , 0];
% P{7}  = [L-ls+w/2, 0];
% P{8}  = [L       , 0];
% P{9}  = [L       , H];
% P{10} = [L/2+w/2 , H];
% P{11} = [L/2     , H];
% P{12} = [L/2-w/2 , H];
% P{13} = [0       , H];

P{1} = [0, 0];
P{2} = [L, 0];
P{3} = [L, H];
P{4} = [0, H];

PS{1} = [ls-w/2  , 0];
PS{2} = [ls      , 0];
PS{3} = [ls+w/2  , 0];
PS{4} = [L-ls-w/2, 0];
PS{5} = [L-ls    , 0];
PS{6} = [L-ls+w/2, 0];
PS{7} = [L/2+w/2 , H];
PS{8} = [L/2     , H];
PS{9} = [L/2-w/2 , H];

C = LINE([b,0],[b,a]);

if refinecrack
    clCrack = clC;
else
    clCrack = [clD clC];
end

numpoints = 1:15;
numlinecrack = 1;
numlines = 2:15;
numlineloop = 1;
numpointscrack = numpoints([2,1]);
numsurface = 1;
G = gmshfile(C,clCrack,numpointscrack,numlinecrack);
G = createpoints(G,P,clD,numpoints([12,6,7,11]));
G = createpoints(G,PS,clS,numpoints([13:15,3:5,8:10]));
G = createcontour(G,numpoints(2:15),numlines,numlineloop);
G = createplanesurface(G,numlineloop,numsurface);
G = embedcurveinsurface(G,numlinecrack,numsurface);

if recombine
    G = recombinesurface(G,numsurface);
end

if ~noduplicate
    openboundaryphysicalgroup = 1;
    physicalgroup = 1;
    G = createphysicalpoint(G,numpointscrack(1),openboundaryphysicalgroup);
    G = createphysicalcurve(G,numlinecrack,physicalgroup);
end

numphysicalsurface = 1;
G = createphysicalsurface(G,numsurface,numphysicalsurface);

% Box field
if ~isempty(Box) && isstruct(Box)
    if isfield(Box,'VIn')
        VIn = Box.VIn;
    else
        VIn = clC;
    end
    if isfield(Box,'VOut')
        VOut = Box.VOut;
    else
        VOut = clD;
    end
    XMin = Box.XMin;
    XMax = Box.XMax;
    YMin = Box.YMin;
    YMax = Box.YMax;
    if indim==3 || isfield(Box,'ZMin')
        ZMin = Box.ZMin;
    else
        ZMin = 0;
    end
    if indim==3 || isfield(Box,'ZMax')
        ZMax = Box.ZMax;
    else
        ZMax = 0;
    end
    if isfield(Box,'Thickness')
        Thickness = Box.Thickness;
    else
        Thickness = 0;
    end
    G = createboxfield(G,VIn,VOut,XMin,XMax,YMin,YMax,ZMin,ZMax,Thickness);
    G = setbgfield(G);
end

if nargin>=11 && ischar(filename)
    G = setfile(G,filename);
end

n = max(nargout,1);
varargout = cell(1,n);
[varargout{:}] = gmsh2femobject(indim,G,dim:-1:dim-n+1,varargin{:});

if ~noduplicate
    G = createcrack(G,dim-1,physicalgroup,openboundaryphysicalgroup);
    G = remesh(G,dim,varargin{:});
    
    [varargout{:}] = gmsh2femobject(indim,G,dim:-1:dim-n+1,varargin{:});
end
