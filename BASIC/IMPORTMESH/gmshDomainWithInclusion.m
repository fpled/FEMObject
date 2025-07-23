function varargout = gmshDomainWithInclusion(D,I,clD,clI,filename,indim,varargin)
% function varargout = gmshDomainWithInclusion(D,I,clD,clI,filename,indim)
% D : DOMAIN or QUADRANGLE
% I : DOMAIN or CIRCLE or ELLIPSE or QUADRANGLE or LINE or POINT
% clD, clI : characteristic lengths
% filename : file name (optional)
% indim : space dimension (optional, getindim(D) by default)

isextruded = ischarin('extrude',varargin);
varargin = delonlycharin('extrude',varargin);

if nargin<6 || isempty(indim)
    indim = getindim(D);
end
if nargin<4 || isempty(clI)
    clI = clD;
end

if ~iscell(I)
    I = {I};
end
if isscalar(clI)
    clI = repmat(clI,1,length(I));
end

dim = getdim(D);

if dim==2
    numpoints = 1:4;
    numlines = 1:4;
    numlineloop = 1;
    G = gmshfile(D,clD,numpoints,numlines,numlineloop);
elseif dim==3
    numpoints = 1:8;
    numlines = 1:12;
    numlineloop = 1:6;
    numsurface = 1:6;
    numsurfaceloop = 1;
    G = gmshfile(D,clD,numpoints,numlines,numlineloop,numsurface,numsurfaceloop);
end
if nargin>=5 && ischar(filename)
    G = setfile(G,filename);
end

numembeddedpoints = [];
numembeddedlines = [];
if dim==2
    numcurves = numlines;
    numsurface = 1;
    numphysicalsurface = 1;
elseif dim==3
    numembeddedsurfaces = [];
    numsurfaces = numsurface;
    numvolume = 1;
    numphysicalvolume = 1;
end
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
            if getdim(I{j})==2
                numpoints = numpoints(end)+(1:4);
                numlines = numlines(end)+(1:4);
                numlineloop = numlineloop(end)+1;
                numsurface = numsurface(end)+1;
                if dim==2
                    GI = gmshfile(I{j},clI(j),numpoints,numlines,numlineloop,numsurface);
                elseif dim==3
                    GI = gmshfile(I{j},clI(j),numpoints,numlines,numlineloop,numsurface);
                    numembeddedsurfaces = [numembeddedsurfaces,numsurface];
                end
            elseif getdim(I{j})==3
                numpoints = numpoints(end)+(1:8);
                numlines = numlines(end)+(1:12);
                numlineloop = numlineloop(end)+(1:6);
                numsurface = numsurface(end)+(1:6);
                numsurfaceloop = numsurfaceloop(end)+1;
                numvolume = numvolume(end)+1;
                GI = gmshfile(I{j},clI(j),numpoints,numlines,numlineloop,numsurface,numsurfaceloop,numvolume);
            end
        elseif isa(I{j},'CIRCLE') || isa(I{j},'ELLIPSE')
            numpoints = numpoints(end)+(1:5);
            numlines = numlines(end)+(1:4);
            numlineloop = numlineloop(end)+1;
            numsurface = numsurface(end)+1;
            if dim==2
                GI = gmshfile(I{j},clI(j),numpoints(1),numpoints(2:end),numlines,numlineloop,numsurface);
            elseif dim==3
                GI = gmshfile(I{j},clI(j),numpoints(1),numpoints(2:end),numlines,numlineloop,numsurface);
                numembeddedsurfaces = [numembeddedsurfaces,numsurface];
            end
        elseif isa(I{j},'SPHERE') || isa(I{j},'ELLIPSOID')
            numpoints = numpoints(end)+(1:7);
            numlines = numlines(end)+(1:12);
            numlineloop = numlineloop(end)+(1:8);
            numsurface = numsurface(end)+(1:8);
            numsurfaceloop = numsurfaceloop(end)+1;
            numvolume = numvolume(end)+1;
            GI = gmshfile(I{j},clI(j),numpoints(1),numpoints(2:end),numlines,numlineloop,numsurface,numsurfaceloop,numvolume);
        elseif isa(I{j},'CYLINDER')
            P = getvertices(I{j});
            tol = getfemobjectoptions('tolerancepoint');
            angle = getangle(I{j});
            n = numel(P)/2; % number of points at base/top
            numpoints = numpoints(end)+(1:(2+2*n)); % 1 base center + 1 top center + n base points + n top points
            if abs(angle - 2*pi) < tol
                % Full cylinder: closed contour
                numlines = numlines(end)+(1:3*n); % n base arcs + n top arcs + n verticals
                numlineloop = numlineloop(end)+(1:(n+2)); % 1 base circle + 1 top circle + n lateral faces
                numsurface = numsurface(end)+(1:(n+2)); % same as numlineloop
            else
                % Partial cylinder: open arc
                numlines = numlines(end)+(1:(3*n+3)); % (n-1) base arcs + 2 base radials + (n-1) top arcs + 2 top radials + n verticals + 1 center vertical
                numlineloop = numlineloop(end)+(1:(n+3)); % 1 base circle + 1 top circle + (n-1) lateral faces +  2 radial faces
                numsurface = numsurface(end)+(1:(n+3)); % same as numlineloop
            end
            numsurfaceloop = numsurfaceloop(end)+1;
            numvolume = numvolume(end)+1;
            GI = gmshfile(I{j},clI(j),numpoints(1:2),numpoints(3:end),numlines,numlineloop,numsurface,numsurfaceloop,numvolume);
        end
        if dim==2 && getdim(I{j})==2
            numphysicalsurface = numphysicalsurface(end)+1;
            GI = createphysicalsurface(GI,numsurface,numphysicalsurface);
            numcurves = [numcurves,-numlines];
        elseif dim==3 && getdim(I{j})==3
            numphysicalvolume = numphysicalvolume(end)+1;
            GI = createphysicalvolume(GI,numvolume,numphysicalvolume);
            numsurfaces = [numsurfaces,-numsurface];
        end
    end
    G = G+GI;
end

if dim==2
    numlineloop = numlineloop(end)+1;
    numsurface = 1;
    G = createcurveloop(G,numcurves,numlineloop);
    G = createplanesurface(G,numlineloop,numsurface);
elseif dim==3
    if isextruded
        P = getvertices(D);
        vect = P{5}-P{1};
        [G,out] = extrude(G,vect,'Surface',1,varargin{:});
        numvolume = [out,'[1]'];
    else
        numsurfaceloop = numsurfaceloop(end)+1;
        numvolume = 1;
        G = createsurfaceloop(G,numsurfaces,numsurfaceloop);
        G = createvolume(G,numsurfaceloop,numvolume);
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

if dim==3 && ~isempty(numembeddedsurfaces)
    G = embedsurfacesinvolume(G,numembeddedsurfaces,numvolume);
end

if ischarin('recombine',varargin)
    G = recombinesurface(G);
end

if dim==2
    numphysicalsurface = 1;
    G = createphysicalsurface(G,numsurface,numphysicalsurface);
elseif dim==3
    numphysicalvolume = 1;
    G = createphysicalvolume(G,numvolume,numphysicalvolume);
end

varargin = delonlycharin('recombine',varargin);

% Box field
B = getcharin('Box',varargin,[]);
if ~isempty(B) && isstruct(B)
    if isfield(B,'VIn')
        VIn = B.VIn;
    else
        VIn = min(clI);
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

n = max(nargout,1);
varargout = cell(1,n);
[varargout{:}] = gmsh2femobject(indim,G,dim:-1:dim-n+1,varargin{:});
