function varargout = gmshCircleWithHole(C,H,clC,clH,filename,indim,varargin)
% function varargout = gmshCircleWithHole(C,H,clC,clH,filename,indim)
% C : CIRCLE
% H : DOMAIN or CIRCLE or ELLIPSE or QUADRANGLE or LIGNE or POINT
% clC, clH : characteristic lengths
% filename : file name (optional)
% indim : space dimension (optional, getindim(D) by default)

noduplicate = ischarin('noduplicate',varargin);
varargin = delonlycharin('noduplicate',varargin);

if nargin<6 || isempty(indim)
    indim = getindim(C);
end
if nargin<4 || isempty(clH)
    clH = clC;
end

if ~iscell(H)
    H = {H};
end
if isscalar(clH)
    clH = repmat(clH,1,length(H));
end

G = gmshfile(C,clC,1,2:5,1:4,5);
if nargin>=5 && ischar(filename)
    G = setfile(G,filename);
end

numpoints = 5+(1:5);
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
if ~isempty(numberembeddedpoints)
    G = embedpointsinsurface(G,numberembeddedpoints,1);
end
if ~isempty(numberembeddedlines)
    G = embedcurvesinsurface(G,numberembeddedlines,1);
    if ~noduplicate
        physicalgroup = 1;
        G = createphysicalpoint(G,numberpointsinembeddedlines,1);
        G = createphysicalcurve(G,numberembeddedlines,physicalgroup);
    end
end
if ischarin('recombine',varargin)
    G = recombinesurface(G,1);
end
G = createphysicalsurface(G,1,1);
varargin = delonlycharin('recombine',varargin);

n=max(nargout,1);
varargout = cell(1,n);
[varargout{:}] = gmsh2femobject(indim,G,getdim(C):-1:getdim(C)-n+1,varargin{:});

if ~noduplicate && ~isempty(numberembeddedlines)
    G = createcrack(G,getdim(C)-1,physicalgroup);
    G = remesh(G,getdim(C),varargin{:});
    G = deleteoptfile(G);
    
    [varargout{:}] = gmsh2femobject(indim,G,getdim(C):-1:getdim(C)-n+1,varargin{:});
end
