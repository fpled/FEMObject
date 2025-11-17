function varargout = gmshAsymmetricPlateWithSingleEdgeCrackThreeHoles(a,b,clD,clC,clH,unit,filename,indim,varargin)
% function varargout = gmshAsymmetricPlateWithSingleEdgeCrackThreeHoles(a,b,clD,clC,clH,unit,filename,indim)
% a : length of the edge crack (or notch)
% b : location of the edge crack (or notch) from the centerline
% clD, clC, clH : characteristic lengths
% unit : unit in m (optional) (1e-3 for mm, 25.4e-3 for inch)
% filename : file name (optional)
% indim : space dimension (optional, 2 by default)

Box         = getcharin('Box',varargin,[]);
noduplicate = ischarin('noduplicate',varargin);
refinecrack = ischarin('refinecrack',varargin);
recombine   = ischarin('recombine',varargin);

varargin = delcharin('Box',varargin);
varargin = delonlycharin({'noduplicate','refinecrack','recombine'},varargin);

if nargin<8 || isempty(indim)
    indim = 2;
end
if nargin<6 || isempty(unit)
    unit = 1e-3;
end
if nargin<5 || isempty(clH)
    clH = clD;
end
if nargin<4 || isempty(clC)
    clC = clD;
end

dim = 2;

L = 10*unit; % half-length
h = 4*unit; % half-height
ls = 9*unit; % location of the support from the centerline
lh = 4*unit; % location of the holes from the centerline
dh = 2*unit; % distance between the holes
ph = 1.25*unit; % location of the top hole from the top
r = 0.25*unit; % radius of the holes

P{1} = [-L,  -h];
P{2} = [-ls, -h];
P{3} = [ ls, -h];
P{4} = [ L , -h];
P{5} = [ L ,  h];
P{6} = [ 0 ,  h];
P{7} = [-L ,  h];

C = LINE([-b,-h],[-b,-h+a]);

H{1} = CIRCLE(-lh, h-ph-2*dh, r);
H{2} = CIRCLE(-lh, h-ph-dh  , r);
H{3} = CIRCLE(-lh, h-ph     , r);

if refinecrack
    clCrack = clC;
else
    clCrack = [clD clC];
end

numpoints = 1:9;
numlinecrack = 1;
numlines = 2:9;
numpointscrack = numpoints([2,1]);
G = gmshfile(C,clCrack,numpointscrack,numlinecrack);
G = createpoints(G,P,clD,numpoints([8,9,3:7]));
G = createcontour(G,numpoints(2:9),numlines);

numlineloop = 0;
numcurves = numlines;
for j=1:length(H)
    numpoints = numpoints(end)+(1:5);
    numlines = numlines(end)+(1:4);
    numlineloop = numlineloop(end)+1;
    GH = gmshfile(H{j},clH,numpoints(1),numpoints(2:end),numlines,numlineloop);
    G = G+GH;
    numcurves = [numcurves,-numlines];
end

numlineloop = numlineloop(end)+1;
numsurface = 1;
G = createcurveloop(G,numcurves,numlineloop);
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

if nargin>=7 && ischar(filename)
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
