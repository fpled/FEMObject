function varargout = gmshDomainWithHole(D,H,clD,clH,filename,indim,varargin)
% function varargout = gmshDomainWithHole(D,H,clD,clH,filename,indim)
% D : DOMAIN or QUADRANGLE
% H : DOMAIN or CIRCLE or ELLIPSE or QUADRANGLE or LINE or POINT
% clD, clH : characteristic lengths
% filename : file name (optional)
% indim : space dimension (optional, getindim(D) by default)

noduplicate = ischarin('noduplicate',varargin);
isextruded = ischarin('extrude',varargin);
varargin = delonlycharin({'noduplicate','extrude'},varargin);

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
numpointsinembeddedlines = [];
numembeddedlines = [];
if dim==2
    numcurves = numlines;
elseif dim==3
    numpointsinembeddedsurfaces = [];
    numlinesinembeddedsurfaces = [];
    numembeddedsurfaces = [];
    numsurfaces = numsurface;
end
for j=1:length(H)
    if isa(H{j},'POINT')
        numpoints = numpoints(end)+1;
        GH = gmshfile(H{j},clH(j),numpoints);
        numembeddedpoints = [numembeddedpoints,numpoints];
    elseif isa(H{j},'LINE')
        numpoints = numpoints(end)+(1:2);
        numlines = numlines(end)+1;
        GH = gmshfile(H{j},clH(j),numpoints,numlines);
        numpointsinembeddedlines = [numpointsinembeddedlines,numpoints];
        numembeddedlines = [numembeddedlines,numlines];
    else
        if isa(H{j},'DOMAIN') || isa(H{j},'QUADRANGLE')
            if getdim(H{j})==2
                numpoints = numpoints(end)+(1:4);
                numlines = numlines(end)+(1:4);
                numlineloop = numlineloop(end)+1;
                if dim==2
                    GH = gmshfile(H{j},clH(j),numpoints,numlines,numlineloop);
                elseif dim==3
                    numsurface = numsurface(end)+1;
                    GH = gmshfile(H{j},clH(j),numpoints,numlines,numlineloop,numsurface);
                    numpointsinembeddedsurfaces = [numpointsinembeddedsurfaces,numpoints];
                    numlinesinembeddedsurfaces = [numlinesinembeddedsurfaces,numlines];
                    numembeddedsurfaces = [numembeddedsurfaces,numsurface];
                end
            elseif getdim(H{j})==3
                numpoints = numpoints(end)+(1:8);
                numlines = numlines(end)+(1:12);
                numlineloop = numlineloop(end)+(1:6);
                numsurface = numsurface(end)+(1:6);
                numsurfaceloop = numsurfaceloop(end)+1;
                GH = gmshfile(H{j},clH(j),numpoints,numlines,numlineloop,numsurface,numsurfaceloop);
            end
        elseif isa(H{j},'CIRCLE') || isa(H{j},'ELLIPSE')
            numpoints = numpoints(end)+(1:5);
            numlines = numlines(end)+(1:4);
            numlineloop = numlineloop(end)+1;
            if dim==2
                GH = gmshfile(H{j},clH(j),numpoints(1),numpoints(2:end),numlines,numlineloop);
            elseif dim==3
                numsurface = numsurface(end)+1;
                GH = gmshfile(H{j},clH(j),numpoints(1),numpoints(2:end),numlines,numlineloop,numsurface);
                numpointsinembeddedsurfaces = [numpointsinembeddedsurfaces,numpoints];
                numlinesinembeddedsurfaces = [numlinesinembeddedsurfaces,numlines];
                numembeddedsurfaces = [numembeddedsurfaces,numsurface];
            end
        elseif isa(H{j},'SPHERE') || isa(H{j},'ELLIPSOID')
            numpoints = numpoints(end)+(1:7);
            numlines = numlines(end)+(1:12);
            numlineloop = numlineloop(end)+(1:8);
            numsurface = numsurface(end)+(1:8);
            numsurfaceloop = numsurfaceloop(end)+1;
            GH = gmshfile(H{j},clH(j),numpoints(1),numpoints(2:end),numlines,numlineloop,numsurface,numsurfaceloop);
        elseif isa(H{j},'CYLINDER')
            P = getvertices(H{j});
            tol = getfemobjectoptions('tolerancepoint');
            angle = getangle(H{j});
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
            GH = gmshfile(H{j},clH(j),numpoints(1:2),numpoints(3:end),numlines,numlineloop,numsurface,numsurfaceloop);
        end
        if dim==2 && getdim(H{j})==2
            numcurves = [numcurves,-numlines];
        elseif dim==3 && getdim(H{j})==3
            numsurfaces = [numsurfaces,-numsurface];
        end
    end
    G = G+GH;
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
        if ~noduplicate
            physicalgroup = 1;
            % G = createphysicalpoint(G,numpointsinembeddedlines,1);
            G = createphysicalcurve(G,numembeddedlines,physicalgroup);
        end
    elseif dim==3
        G = embedcurvesinvolume(G,numembeddedlines,numvolume);
    end
end

if dim==3 && ~isempty(numembeddedsurfaces)
    G = embedsurfacesinvolume(G,numembeddedsurfaces,numvolume);
    if ~noduplicate
        physicalgroup = 1;
        % G = createphysicalpoint(G,numpointsinembeddedsurfaces,1);
        % G = createphysicalcurve(G,numlinesinembeddedsurfaces,1);
        G = createphysicalsurface(G,numembeddedsurfaces,physicalgroup);
    end
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

n = max(nargout,1);
varargout = cell(1,n);
[varargout{:}] = gmsh2femobject(indim,G,dim:-1:dim-n+1,varargin{:});

if ~noduplicate && ( (dim==2 && ~isempty(numembeddedlines)) || (dim==3 && ~isempty(numembeddedsurfaces)) )
    G = createcrack(G,dim-1,physicalgroup);
    G = remesh(G,dim,varargin{:});
    G = deleteoptfile(G);
    
    [varargout{:}] = gmsh2femobject(indim,G,dim:-1:dim-n+1,varargin{:});
end
