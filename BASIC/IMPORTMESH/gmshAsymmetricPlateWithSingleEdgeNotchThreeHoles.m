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

H{1} = CIRCLE(-lh,h-ph-2*dh,r);
H{2} = CIRCLE(-lh,h-ph-dh,r);
H{3} = CIRCLE(-lh,h-ph,r);

G = GMSHFILE();
if nargin>=8 && ischar(filename)
    G = setfile(G,filename);
end

if ischarin('r',varargin) || ischarin('rect',varargin) || ischarin('rectangular',varargin)
    % rectangular notch
    PC{1} = [-b-c/2,-h];
    PC{2} = [-b-c/2,-h+a];
    PC{3} = [-b+c/2,-h+a];
    PC{4} = [-b+c/2,-h];
    G = createpoints(G,P,clD,[1:2,7:11]);
    if ischarin('refinecrack',varargin)
        clcrack = clC;
    else
        clcrack = [clD clC clC clD];
    end
    G = createpoints(G,PC,clcrack,3:6);
    G = createlines(G,[[1 2];[2 3];[3 4];[4 5];[5 6];...
        [6 7];[7 8];[8 9];[9 10];[10 11];[11 1]],1:11);
elseif ischarin('v',varargin) || ischarin('V',varargin)|| ischarin('triangular',varargin)
    % V (triangular) notch
    PC{1} = [-b-c/2,-h];
    PC{2} = [-b,-h+a];
    PC{3} = [-b+c/2,-h];
    G = createpoints(G,P,clD,[1:2,6:10]);
    if ischarin('refinecrack',varargin)
        clcrack = clC;
    else
        clcrack = [clD clC clD];
    end
    G = createpoints(G,PC,clcrack,3:5);
    G = createlines(G,[[1 2];[2 3];[3 4];[4 5];[5 6];...
        [6 7];[7 8];[8 9];[9 10];[10 1]],1:10);
else%if ischarin('c',varargin) || ischarin('circ',varargin) || ischarin('circular',varargin)
    % circular notch
    PC{1} = [-b-c/2,-h];
    PC{2} = [-b-c/2,-h+a-c/2];
    PC{3} = [-b,-h+a];
    PC{4} = [-b+c/2,-h+a-c/2];
    PC{5} = [-b,-h+a-c/2];
    PC{6} = [-b+c/2,-h];
    G = createpoints(G,P,clD,[1:2,9:13]);
    if ischarin('refinecrack',varargin)
        clcrack = clC;
    else
        clcrack = [clD clC clC clC clC clD];
    end
    G = createpoints(G,PC,clcrack,3:8);
    G = createcircle(G,7,4:5,4);
    G = createcircle(G,7,5:6,5);
    G = createlines(G,[[1 2];[2 3];[3 4];...
        [6 8];[8 9];[9 10];[10 11];[11 12];[12 13];[13 1]],[1:3,6:12]);
end

if ischarin('r',varargin) || ischarin('rect',varargin) || ischarin('rectangular',varargin)
    % rectangular notch
    numpoints = 12:16;
    numlines = 12:16;
    numlineloop = 1:11;
elseif ischarin('v',varargin) || ischarin('V',varargin)|| ischarin('triangular',varargin)
    % V (triangular) notch
    numpoints = 11:15;
    numlines = 11:15;
    numlineloop = 1:10;
else%if ischarin('c',varargin) || ischarin('circ',varargin) || ischarin('circular',varargin)
    % circular notch
    numpoints = 14:18;
    numlines = 13:17;
    numlineloop = 1:12;
end
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
if ischarin('r',varargin) || ischarin('rect',varargin) || ischarin('rectangular',varargin)
    % rectangular notch
    numlinecrack = 3:5;
elseif ischarin('v',varargin) || ischarin('V',varargin)|| ischarin('triangular',varargin)
    % V (triangular) notch
    numlinecrack = 3:4;
else%if ischarin('c',varargin) || ischarin('circ',varargin) || ischarin('circular',varargin)
    % circular notch
    numlinecrack = 3:6;
end
G = createphysicalcurve(G,numlinecrack,1);
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

n=max(nargout,1);
varargout = cell(1,n);
[varargout{:}] = gmsh2femobject(indim,G,2:-1:2-n+1,varargin{:});
