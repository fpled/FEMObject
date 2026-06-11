function varargout = gmshThreePointBendingWithSingleEdgeNotch(L,H,ls,w,a,b,c,t,clD,clC,filename,indim,varargin)
% function varargout = gmshThreePointBendingWithSingleEdgeNotch(L,H,ls,w,a,b,c,t,clD,clC,filename,indim,notchtype)
% L : length
% H : height
% ls : location of the supports from the outer edges
% w : flat punch width
% a : length of the edge crack (or notch)
% b : location of the edge crack (or notch) from the left outer edge
% c : width of the edge crack (or notch)
% t : thickness
% clD, clC : characteristic lengths
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

if nargin<12 || isempty(indim)
    indim = 2;
end
if nargin<10 || isempty(clC)
    clC = clD;
end
if nargin<8 || isempty(t)
    t = 1;
end

dim = 2;

P{1}  = [0       , 0];
P{2}  = [ls-w/2  , 0];
P{3}  = [ls      , 0];
P{4}  = [ls+w/2  , 0];
P{5}  = [L-ls-w/2, 0];
P{6}  = [L-ls    , 0];
P{7}  = [L-ls+w/2, 0];
P{8}  = [L       , 0];
P{9}  = [L       , H];
P{10} = [L/2+w/2 , H];
P{11} = [L/2     , H];
P{12} = [L/2-w/2 , H];
P{13} = [0       , H];

G = GMSHFILE();
if nargin>=11 && ischar(filename)
    G = setfile(G,filename);
end

if isrect
    % rectangular notch
    PC{1} = [b-c/2, 0];
    PC{2} = [b-c/2, a];
    PC{3} = [b+c/2, a];
    PC{4} = [b+c/2, 0];
    
    if refinecrack
        clCrack = clC;
    else
        clCrack = [clD clC clC clD];
    end
    
    numpoints = 1:17;
    numlines  = 1:17;
    G = createpoints(G,P,clD,numpoints([12:15,3:13]));
    G = createpoints(G,PC,clCrack,numpoints([16:17,1:2]));
    n = length(numpoints);
    seg = [1:n; 2:n,1]';
    G = createlines(G,numpoints(seg),numlines);
elseif istri
    % V (triangular) notch
    PC{1} = [b-c/2, 0];
    PC{2} = [b    , a];
    PC{3} = [b+c/2, 0];
    
    if refinecrack
        clCrack = clC;
    else
        clCrack = [clD clC clD];
    end
    
    numpoints = 1:16;
    numlines  = 1:16;
    G = createpoints(G,P,clD,numpoints([12:15,3:11]));
    G = createpoints(G,PC,clCrack,numpoints([16,1:2]));
    n = length(numpoints);
    seg = [1:n; 2:n,1]';
    G = createlines(G,numpoints(seg),numlines);
else%if iscirc
    % circular notch
    PC{1} = [b-c/2, 0    ];
    PC{2} = [b-c/2, a-c/2];
    PC{3} = [b    , a    ];
    PC{4} = [b+c/2, a-c/2];
    PC{5} = [b    , a-c/2];
    PC{6} = [b+c/2, 0    ];
    
    if refinecrack
        clCrack = clC;
    else
        clCrack = [clD clC clC clC clC clD];
    end
    
    numpoints = 1:19;
    numlines  = 1:18;
    G = createpoints(G,P,clD,numpoints([12:15,3:11]));
    G = createpoints(G,PC,clCrack,numpoints([16:18,1,19,2]));
    G = createcirclearc(G,numpoints(19),numpoints([17,18]),numlines(1));
    G = createcirclearc(G,numpoints(19),numpoints([18,1]),numlines(2));
    G = createlines(G,numpoints([1:16;2:17]'),numlines(3:18));
end

numlineloop = 1;
numsurface = 1;
G = createcurveloop(G,numlines,numlineloop);
G = createplanesurface(G,numlineloop,numsurface);

if recombine
    G = recombinesurface(G,numsurface);
end

if isrect
    % rectangular notch
    numlinecrack = [16:17,1];
elseif istri
    % V (triangular) notch
    numlinecrack = [16,1];
else%if iscirc
    % circular notch
    numlinecrack = [18,1:3];
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
