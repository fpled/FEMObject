function varargout = gmshAsymmetricPlateWithSingleEdgeNotchThreeHoles(a,b,c,clD,clC,clH,unit,filename,indim,varargin)
% function varargout = gmshAsymmetricPlateWithSingleEdgeNotchThreeHoles(a,b,c,clD,clC,clH,unit,filename,indim,notchtype)
% a : length of the edge crack (or notch)
% b : location of the edge crack (or notch) from the centerline
% c : width of the edge crack (or notch)
% clD, clC, clH : characteristic lengths
% unit : unit in m (optional) (1e-3 for mm, 25.4e-3 for inch)
% filename : file name (optional)
% indim : space dimension (optional, 2 by default)
% notchtype : 'c', 'circ' or 'circular' for circular notch (by default)
%             'r', 'rect' or 'rectangular' for rectangular notch
%             'v', 'V' or 'triangular' for V notch

Box         = getcharin('Box',varargin,[]);
refinecrack = ischarin('refinecrack',varargin);
recombine   = ischarin('recombine',varargin);
isrect      = any(ischarin({'r','rect','rectangular'},varargin));
istri       = any(ischarin({'v','V','triangular'},varargin));
iscirc      = any(ischarin({'c','circ','circular'},varargin));

varargin = delcharin('Box',varargin);
varargin = delonlycharin({'refinecrack','recombine', ...
                          'r','rect','rectangular', ...
                          'v','V','triangular', ...
                          'c','circ','circular'},varargin);

if nargin<9 || isempty(indim)
    indim = 2;
end
if nargin<7 || isempty(unit)
    unit = 1;
end
if nargin<6 || isempty(clH)
    clH = clD;
end
if nargin<5 || isempty(clC)
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

P{1} = [-L , -h];
P{2} = [-ls, -h];
P{3} = [ ls, -h];
P{4} = [ L , -h];
P{5} = [ L ,  h];
P{6} = [ 0 ,  h];
P{7} = [-L ,  h];

H{1} = CIRCLE(-lh, h-ph-2*dh, r);
H{2} = CIRCLE(-lh, h-ph-dh  , r);
H{3} = CIRCLE(-lh, h-ph     , r);

G = GMSHFILE();
if nargin>=8 && ischar(filename)
    G = setfile(G,filename);
end

if isrect
    % rectangular notch
    PC{1} = [-b-c/2, -h  ];
    PC{2} = [-b-c/2, -h+a];
    PC{3} = [-b+c/2, -h+a];
    PC{4} = [-b+c/2, -h  ];
    
    if refinecrack
        clCrack = clC;
    else
        clCrack = [clD clC clC clD];
    end
    
    numpoints = 1:11;
    numlines  = 1:11;
    G = createpoints(G,P,clD,numpoints([1:2,7:11]));
    G = createpoints(G,PC,clCrack,numpoints(3:6));
    n = length(numpoints);
    seg = [1:n; 2:n,1]';
    G = createlines(G,numpoints(seg),numlines);
elseif istri
    % V (triangular) notch
    PC{1} = [-b-c/2, -h  ];
    PC{2} = [-b    , -h+a];
    PC{3} = [-b+c/2, -h  ];
    
    if refinecrack
        clCrack = clC;
    else
        clCrack = [clD clC clD];
    end
    
    numpoints = 1:10;
    numlines  = 1:10;
    G = createpoints(G,P,clD,numpoints([1:2,6:10]));
    G = createpoints(G,PC,clCrack,numpoints(3:5));
    n = length(numpoints);
    seg = [1:n; 2:n,1]';
    G = createlines(G,numpoints(seg),numlines);
else%if iscirc
    % circular notch
    PC{1} = [-b-c/2, -h      ];
    PC{2} = [-b-c/2, -h+a-c/2];
    PC{3} = [-b    , -h+a    ];
    PC{4} = [-b+c/2, -h+a-c/2];
    PC{5} = [-b    , -h+a-c/2];
    PC{6} = [-b+c/2, -h      ];
    
    if refinecrack
        clCrack = clC;
    else
        clCrack = [clD clC clC clC clC clD];
    end
    
    numpoints = 1:13;
    numlines  = 1:12;
    G = createpoints(G,P,clD,numpoints([1:2,9:13]));
    G = createpoints(G,PC,clCrack,numpoints(3:8));
    G = createcirclearc(G,numpoints(7),numpoints(4:5),numlines(4));
    G = createcirclearc(G,numpoints(7),numpoints(5:6),numlines(5));
    seg = [1:3, 6, 8:13; 2:4, 8:13, 1]';
    G = createlines(G,numpoints(seg),numlines([1:3,6:12]));
end

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

if recombine
    G = recombinesurface(G,numsurface);
end

if isrect
    % rectangular notch
    numlinecrack = 3:5;
elseif istri
    % V (triangular) notch
    numlinecrack = 3:4;
else%if iscirc
    % circular notch
    numlinecrack = 3:6;
end
% numphysicalcurve = 1;
numphysicalsurface = 1;
% G = createphysicalcurve(G,numlinecrack,numphysicalcurve);
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

n = max(nargout,1);
varargout = cell(1,n);
[varargout{:}] = gmsh2femobject(indim,G,dim:-1:dim-n+1,varargin{:});
