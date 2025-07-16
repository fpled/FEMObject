function varargout = gmshDomainWithHole(D,H,clD,clH,filename,indim,varargin)
% function varargout = gmshDomainWithHole(D,H,clD,clH,filename,indim)
% D : DOMAIN or QUADRANGLE
% H : DOMAIN or CIRCLE or ELLIPSE or QUADRANGLE or LIGNE or POINT
% clD, clH : characteristic lengths
% filename : file name (optional)
% indim : space dimension (optional, getindim(D) by default)

noduplicate = ischarin('noduplicate',varargin);
varargin = delonlycharin('noduplicate',varargin);

if nargin<6 || isempty(indim)
    indim = getindim(D);
end
if nargin<4 || isempty(clH)
    clH = clD;
end

if ~iscell(H)
    H = {H};
end
if isscalar(clH)
    clH = repmat(clH,1,length(H));
end

G = gmshfile(D,clD,1:4,1:4,5);
if nargin>=5 && ischar(filename)
    G = setfile(G,filename);
end

numpoints = 4+(1:5);
numlines = 5+(1:5);
numlineloop = 1:4;
numberembeddedpoints = [];
numberpointsinembeddedlines = [];
numberembeddedlines = [];
for j=1:length(H)
    if isa(H{j},'POINT')
        GH = gmshfile(H{j},clH(j),numpoints(1));
        numberembeddedpoints = [numberembeddedpoints,numpoints(1)];
    elseif isa(H{j},'LIGNE')
        GH = gmshfile(H{j},clH(j),numpoints(1:2),numlines(1));
        numberpointsinembeddedlines = [numberpointsinembeddedlines,numpoints(1:2)];
        numberembeddedlines = [numberembeddedlines,numlines(1)];
    else
        if isa(H{j},'DOMAIN') || isa(H{j},'QUADRANGLE')
            GH = gmshfile(H{j},clH(j),numpoints(1:end-1),numlines(1:end-1),numlines(end));
        elseif isa(H{j},'CIRCLE') || isa(H{j},'ELLIPSE')
            GH = gmshfile(H{j},clH(j),numpoints(1),numpoints(2:end),numlines(1:end-1),numlines(end));
        end
        numlineloop = [numlineloop,-numlines(1:end-1)];
    end
    G = G+GH;
    numpoints = numpoints+5;
    numlines = numlines+5;
end
G = createcurveloop(G,numlineloop,numlines(end));
G = createplanesurface(G,numlines(end),1);
if getdim(D)==3
    P = getvertices(D);
    vect = P{5}-P{1};
    G = extrude(G,vect,'Surface',1,varargin{:});
end
if ischarin('recombine',varargin)
    if getdim(D)==2
        G = recombinesurface(G,1);
    elseif getdim(D)==3
        G = recombinesurface(G);
    end
end
if ~isempty(numberembeddedpoints)
    if getdim(D)==2
        G = embedpointsinsurface(G,numberembeddedpoints,1);
    elseif getdim(D)==3
        G = embedpointsinvolume(G,numberembeddedpoints,1);
    end
end
if ~isempty(numberembeddedlines)
    if getdim(D)==2
        G = embedcurvesinsurface(G,numberembeddedlines,1);
    elseif getdim(D)==3
        G = embedcurvesinvolume(G,numberembeddedlines,1);
    end
    if ~noduplicate
        physicalgroup = 1;
        G = createphysicalpoint(G,numberpointsinembeddedlines,1);
        G = createphysicalcurve(G,numberembeddedlines,physicalgroup);
    end
end
G = createphysicalsurface(G,1,1);
if getdim(D)==3
    G = createphysicalvolume(G,1,1);
end
varargin = delonlycharin('recombine',varargin);

% Box field
B = getcharin('Box',varargin,[]);
if ~isempty(B) && isstruct(B)
    if isfield(B,'VIn')
        VIn = B.VIn;
    else
        VIn = min(clH);
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
[varargout{:}] = gmsh2femobject(indim,G,getdim(D):-1:getdim(D)-n+1,varargin{:});

if ~noduplicate && ~isempty(numberembeddedlines)
    G = createcrack(G,getdim(D)-1,physicalgroup);
    G = remesh(G,getdim(D),varargin{:});
    G = deleteoptfile(G);
    
    [varargout{:}] = gmsh2femobject(indim,G,getdim(D):-1:getdim(D)-n+1,varargin{:});
end
