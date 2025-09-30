function varargout = gmshCircleWithInclusion(C,I,clC,clI,filename,indim,varargin)
% function varargout = gmshCircleWithInclusion(C,I,clC,clI,filename,indim)
% C : CIRCLE or CYLINDER
% I : DOMAIN or CIRCLE or ELLIPSE or QUADRANGLE or LINE or POINT
% clC, clI : characteristic lengths
% filename : file name (optional)
% indim : space dimension (optional, getindim(C) by default)

Box = getcharin('Box',varargin,[]);
varargin = delcharin('Box',varargin);

if nargin<6 || isempty(indim)
    indim = getindim(C);
end
if nargin<4 || isempty(clI)
    clI = clC;
end

if ~iscell(I)
    I = {I};
end
if isscalar(clI)
    clI = repmat(clI,1,length(I));
end

br = @(tag,k) sprintf('%s[%d]', tag, k); % bracket reference helper

if isa(C,'CIRCLE')
    dim = getdim(C);
elseif isa(C,'CYLINDER')
    dim = getdim(C);
    height = geth(C);
    C = getcircle(C);
end

numpoints = 1:5;
numlines = 1:4;
numlineloop = 1;
G = gmshfile(C,clC,numpoints(1),numpoints(2:end),numlines,numlineloop);
if nargin>=5 && ischar(filename)
    G = setfile(G,filename);
end

numembeddedpoints = [];
numembeddedlines = [];
numcurves = numlines;
numsurface = 1;
numphysicalsurface = 1;
for j=1:length(I)
    if isa(I{j},'POINT')
        numpoints = numpoints(end)+1;
        GI = gmshfile(I{j},clI(j),numpoints);
        numembeddedpoints = [numembeddedpoints,numpoints];
    elseif isa(I{j},'LINE')
        numpoints = numpoints(end)+(1:2);
        numlines = numlines(end)+1;
        GI = gmshfile(I{j},clI(j),numpoints,numlines);
        numembeddedlines = [numembeddedlines,numlines];
    else
        if isa(I{j},'DOMAIN') || isa(I{j},'QUADRANGLE')
            numpoints = numpoints(end)+(1:4);
            numlines = numlines(end)+(1:4);
            numlineloop = numlineloop(end)+1;
            numsurface = numsurface(end)+1;
            GI = gmshfile(I{j},clI(j),numpoints,numlines,numlineloop,numsurface);
        elseif isa(I{j},'CIRCLE') || isa(I{j},'ELLIPSE')
            numpoints = numpoints(end)+(1:5);
            numlines = numlines(end)+(1:4);
            numlineloop = numlineloop(end)+1;
            numsurface = numsurface(end)+1;
            GI = gmshfile(I{j},clI(j),numpoints(1),numpoints(2:end),numlines,numlineloop,numsurface);
        end
        if ischarin('recombine',varargin)
            GI = recombinesurface(GI,numsurface);
        end
        numphysicalsurface = numphysicalsurface(end)+1;
        GI = createphysicalsurface(GI,numsurface,numphysicalsurface);
        numcurves = [numcurves,-numlines];
    end
    G = G+GI;
end

numlineloop = numlineloop(end)+1;
numsurface = 1;
G = createcurveloop(G,numcurves,numlineloop);
G = createplanesurface(G,numlineloop,numsurface);

if dim==3
    normal = getnormal(C);
    axis   = height * normal;
    if ischarin('recombine',varargin) && ~ischarin('Layers',varargin)
        numlayers = max(1, round(height/clC));
        varargin = [varargin, {'Layers',numlayers}];
    end
    [G,tag] = extrude(G,axis,'Surface',numsurface,varargin{:});
    numtopsurface = br(tag,0); % top surface
    numvolume     = br(tag,1); % volume
    if ischarin('recombine',varargin)
        G = recombinesurface(G,{numsurface,numtopsurface});
    end
end

if ~isempty(numembeddedpoints)
    if dim==2
        G = embedpointsinsurface(G,numembeddedpoints,numsurface);
    elseif dim==3
        G = embedpointsinvolume(G,numembeddedpoints,numvolume);
    end
end

if ~isempty(numembeddedlines)
    if dim==2
        G = embedcurvesinsurface(G,numembeddedlines,numsurface);
    elseif dim==3
        G = embedcurvesinvolume(G,numembeddedlines,numvolume);
    end
end

if ischarin('recombine',varargin)
    if dim==2
        G = recombinesurface(G,numsurface);
    elseif dim==3 && ~extrusion
        G = recombinesurface(G);
    end
end

varargin = delonlycharin('recombine',varargin);

if dim==2
    numphysicalsurface = 1;
    G = createphysicalsurface(G,numsurface,numphysicalsurface);
elseif dim==3
    numphysicalvolume = 1;
    G = createphysicalvolume(G,numvolume,numphysicalvolume);
end

% Box field
if ~isempty(Box) && isstruct(Box)
    if isfield(Box,'VIn')
        VIn = Box.VIn;
    else
        VIn = min(clI);
    end
    if isfield(Box,'VOut')
        VOut = Box.VOut;
    else
        VOut = clC;
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
