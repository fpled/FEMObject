function varargout = gmshAsymmetricPlateWithSingleEdgeCrackThreeHoles(a,b,clD,clC,clH,unit,filename,indim,varargin)
% function varargout = gmshAsymmetricPlateWithSingleEdgeCrackThreeHoles(a,b,clD,clC,clH,unit,filename,indim)
% a : length of the edge crack (or notch)
% b : location of the edge crack (or notch) from the centerline
% clD, clC, clH : characteristic lengths
% unit : unit in m (optional) (1e-3 for mm, 25.4e-3 for inch)
% filename : file name (optional)
% indim : space dimension (optional, 2 by default)

noduplicate = ischarin('noduplicate',varargin);
varargin = delonlycharin('noduplicate',varargin);

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

L = 10*unit; % half-length
h = 4*unit; % half-height
ls = 9*unit; % location of the support from the centerline
lh = 4*unit; % location of the holes from the centerline
dh = 2*unit; % distance between the holes
ph = 1.25*unit; % location of the top hole from the top
r = 0.25*unit; % radius of the holes

P{1} = [-L,-h];
P{2} = [-ls,-h];
P{3} = [ls,-h];
P{4} = [L,-h];
P{5} = [L,h];
P{6} = [0,h];
P{7} = [-L,h];

C = LIGNE([-b,-h],[-b,-h+a]);

H{1} = CIRCLE(-lh,h-ph-2*dh,r);
H{2} = CIRCLE(-lh,h-ph-dh,r);
H{3} = CIRCLE(-lh,h-ph,r);

if ischarin('refinecrack',varargin)
    clcrack = clC;
else
    clcrack = [clD clC];
end
G = gmshfile(C,clcrack,[3 9],1);
G = createpoints(G,P,clD,[1:2,4:8]);
G = createcontour(G,1:8,2:9,10);

numpoints = 10:14;
numlines = 11:15;
numlineloop = 2:9;
for j=1:length(H)
    numlineloop = [numlineloop,-numlines(1:end-1)];
    GH = gmshfile(H{j},clH,numpoints(1),numpoints(2:end),numlines(1:end-1),numlines(end));
    G = G+GH;
    numpoints = numpoints+5;
    numlines = numlines+5;
end
G = createcurveloop(G,numlineloop,numlines(end));
G = createplanesurface(G,numlines(end),1);
if ischarin('recombine',varargin)
    G = recombinesurface(G,1);
end
G = embedcurveinsurface(G,1,1);
if ~noduplicate
    openboundaryphysicalgroup = 1;
    physicalgroup = 1;
    G = createphysicalpoint(G,3,openboundaryphysicalgroup);
    G = createphysicalcurve(G,1,physicalgroup);
end
G = createphysicalsurface(G,1,1);

varargin = delonlycharin({'recombine','refinecrack'},varargin);

% Box field
B = getcharin('Box',varargin,[]);
if ~isempty(B) && isstruct(B)
    if isfield(B,'VIn')
        VIn = B.VIn;
    else
        VIn = clC;
    end
    if isfield(B,'VOut')
        VOut = B.VOut;
    else
        VOut = clD;
    end
    XMin = B.XMin;
    XMax = B.XMax;
    YMin = B.YMin;
    YMax = B.YMax;
    if indim==3 || isfield(B,'ZMin')
        ZMin = B.ZMin;
    else
        ZMin = 0;
    end
    if indim==3 || isfield(B,'ZMax')
        ZMax = B.ZMax;
    else
        ZMax = 0;
    end
    if isfield(B,'Thickness')
        Thickness = B.Thickness;
    else
        Thickness = 0;
    end
    G = createboxfield(G,VIn,VOut,XMin,XMax,YMin,YMax,ZMin,ZMax,Thickness);
    G = setbgfield(G);
end

if nargin>=7 && ischar(filename)
    G = setfile(G,filename);
end

n=max(nargout,1);
varargout = cell(1,n);
[varargout{:}] = gmsh2femobject(indim,G,2:-1:2-n+1,varargin{:});

if ~noduplicate
    G = createcrack(G,2-1,physicalgroup,openboundaryphysicalgroup);
    G = remesh(G,2,varargin{:});
    G = deleteoptfile(G);
    
    [varargout{:}] = gmsh2femobject(indim,G,2:-1:2-n+1,varargin{:});
end
