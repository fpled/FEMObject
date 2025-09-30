function varargout = gmshDomainWithHole(D,H,clD,clH,filename,indim,varargin)
% function varargout = gmshDomainWithHole(D,H,clD,clH,filename,indim)
% D : DOMAIN or QUADRANGLE
% H : DOMAIN or CIRCLE or ELLIPSE or QUADRANGLE or LINE or POINT
% clD, clH : characteristic lengths
% filename : file name (optional)
% indim : space dimension (optional, getindim(D) by default)

Box         = getcharin('Box',varargin,[]);
noduplicate = ischarin('noduplicate',varargin);
extrusion   = ischarin('extrude',varargin);

varargin = delcharin('Box',varargin);
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

br = @(tag,k) sprintf('%s[%d]', tag, k); % bracket reference helper

dim = getdim(D);

if dim==2
    numpoints = 1:4;
    numlines = 1:4;
    numlineloop = 1;
    G = gmshfile(D,clD,numpoints,numlines,numlineloop);
elseif dim==3
    if extrusion
        numpoints = 1:4;
        numlines = 1:4;
        numlineloop = 1;
        numsurface = 1;
    else
        numpoints = 1:8;
        numlines = 1:12;
        numlineloop = 1:6;
        numsurface = 1:6;
    end
    numsurfaceloop = 1;
    G = gmshfile(D,clD,numpoints,numlines,numlineloop,numsurface,numsurfaceloop);
end
if nargin>=5 && ischar(filename)
    G = setfile(G,filename);
end

numembeddedpoints = [];
numpointsinembeddedlines = [];
numembeddedlines = [];
if dim==2 || (dim==3 && extrusion)
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
                if dim==2 || (dim==3 && extrusion)
                    GH = gmshfile(H{j},clH(j),numpoints,numlines,numlineloop);
                elseif dim==3 && ~extrusion
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
            if dim==2 || (dim==3 && extrusion)
                GH = gmshfile(H{j},clH(j),numpoints(1),numpoints(2:end),numlines,numlineloop);
            elseif dim==3 && ~extrusion
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
            if isstring(angle), angle = char(angle); end
            if ischar(angle),   angle = str2num(lower(angle)); end
            isfull = abs(angle - 2*pi) < tol;
            n = numel(P)/2; % number of points at base/top
            if extrusion
                numpoints = numpoints(end)+(1:(1+n)); % 1 base center + n base points
                if isfull, numlines = numlines(end)+(1:n); % n base arcs
                else,      numlines = numlines(end)+(1:(n+1)); % (n-1) base arcs + 2 base radials
                end
                numlineloop = numlineloop(end)+1; % 1 base circle
                numsurface = numsurface(end)+1; % same as numlineloop
                center = getc(H{j});
                GH = GMSHFILE();
                GH = createpoint(GH,center,clH(j),numpoints(1));
                GH = createpoints(GH,P(1:n),clH(j),numpoints(2:end));
                if isfull
                    % Full cylinder base: closed circle contour
                    GH = createcirclecontour(GH,numpoints(1),numpoints(2:end),numlines,numlineloop);
                else
                    % Partial cylinder base: open circle arc + two radial lines
                    GH = createcirclearccontour(GH,numpoints(1),numpoints(2:end),numlines,numlineloop);
                end
                H{j} = setdim(H{j},2);
            else
                numpoints = numpoints(end)+(1:(2+2*n)); % 1 base center + 1 top center + n base points + n top points
                if isfull
                    % Full cylinder
                    numlines = numlines(end)+(1:3*n); % n base arcs + n top arcs + n verticals
                    numlineloop = numlineloop(end)+(1:(n+2)); % 1 base circle + 1 top circle + n lateral faces
                    numsurface = numsurface(end)+(1:(n+2)); % same as numlineloop
                else
                    % Partial cylinder
                    numlines = numlines(end)+(1:(3*n+3)); % (n-1) base arcs + 2 base radials + (n-1) top arcs + 2 top radials + n verticals + 1 center vertical
                    numlineloop = numlineloop(end)+(1:(n+3)); % 1 base circle + 1 top circle + (n-1) lateral faces +  2 radial faces
                    numsurface = numsurface(end)+(1:(n+3)); % same as numlineloop
                end
                numsurfaceloop = numsurfaceloop(end)+1;
                GH = gmshfile(H{j},clH(j),numpoints(1:2),numpoints(3:end),numlines,numlineloop,numsurface,numsurfaceloop);
            end
        elseif isa(H{j},'TORUS')
            P = getvertices(H{j});
            tol = getfemobjectoptions('tolerancepoint');
            angle = getangle(H{j});
            if isstring(angle), angle = char(angle); end
            if ischar(angle),   angle = str2num(lower(angle)); end
            isfull = abs(angle - 2*pi) < tol;
            n = numel(P)/4; % number of minor-circle rings
            numpoints = numpoints(end)+(1:(3+5*n)); % 3 major centers + n minor centers + 4*n minor points
            if isfull
                % Full torus
                numlines = numlines(end)+(1:8*n); % 4*n minor arcs + 4*n major arcs
                numlineloop = numlineloop(end)+(1:4*n); % 4*n major circle arcs
                numsurface = numsurface(end)+(1:4*n); % same as numlineloop
            else
                % Partial torus
                numlines = numlines(end)+(1:(8*n-4)); % 4*n minor arcs + 4*(n-1) major arcs
                numlineloop = numlineloop(end)+(1:(4*n-2)); % 4*(n-1) major circle arcs + 2 end minor circles
                numsurface = numsurface(end)+(1:(4*n-2)); % same as numlineloop
            end
            numsurfaceloop = numsurfaceloop(end)+1;
            GH = gmshfile(H{j},clH(j),numpoints(1:3+n),numpoints(3+n+1:end),numlines,numlineloop,numsurface,numsurfaceloop);
        end
        if (dim==2 || (dim==3 && extrusion)) && getdim(H{j})==2
            numcurves = [numcurves,-numlines];
        elseif dim==3 && getdim(H{j})==3
            numsurfaces = [numsurfaces,-numsurface];
        end
    end
    G = G+GH;
end

if dim==2 || (dim==3 && extrusion)
    numlineloop = numlineloop(end)+1;
    numsurface = 1;
    G = createcurveloop(G,numcurves,numlineloop);
    G = createplanesurface(G,numlineloop,numsurface);
elseif dim==3 && ~extrusion
    numsurfaceloop = numsurfaceloop(end)+1;
    numvolume = 1;
    G = createsurfaceloop(G,numsurfaces,numsurfaceloop);
    G = createvolume(G,numsurfaceloop,numvolume);
end

if ~isempty(numembeddedpoints)
    if dim==2 || (dim==3 && extrusion)
        G = embedpointsinsurface(G,numembeddedpoints,numsurface);
    elseif dim==3 && ~extrusion
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
    elseif dim==3 && extrusion
        G = embedcurvesinsurface(G,numembeddedlines,numsurface);
    elseif dim==3 && ~extrusion
        G = embedcurvesinvolume(G,numembeddedlines,numvolume);
    end
end

if dim==3 && ~extrusion && ~isempty(numembeddedsurfaces)
    G = embedsurfacesinvolume(G,numembeddedsurfaces,numvolume);
    if ~noduplicate
        physicalgroup = 1;
        % G = createphysicalpoint(G,numpointsinembeddedsurfaces,1);
        % G = createphysicalcurve(G,numlinesinembeddedsurfaces,1);
        G = createphysicalsurface(G,numembeddedsurfaces,physicalgroup);
    end
end

if dim==3 && extrusion
    % Extrude domain without inclusions
    P = getvertices(D);
    vect = P{5}-P{1};
    if ischarin('recombine',varargin) && ~ischarin('Layers',varargin)
        numlayers = max(1, round(vect(3)/mean([clD,clH(:)'])));
        varargin = [varargin, {'Layers',numlayers}];
    end
    [G,tag] = extrude(G,vect,'Surface',numsurface,varargin{:});
    numtopsurface = br(tag,0); % top surface
    numvolume     = br(tag,1); % volume
    % Extrude embedded lines
    if ~isempty(numembeddedlines)
        N = numel(numembeddedlines);
        % numtoppoints = cell(1,N);
        numtoplines = cell(1,N);
        numsurfaces = cell(1,N);
        for i=1:N
            tag = ['out' num2str(i)];
            G = extrude(G,vect,'Curve',numembeddedlines(i),tag,varargin{:});
            numtoplines{i} = br(tag,0); % top line
            numsurfaces{i} = br(tag,1); % surface
            % numtoppoints{i} = ['p' num2str(i)];
            % G = pointsofcurve(G,numtoplines{i},numtoppoints{i});
        end
        G = embedcurvesinsurface(G,numtoplines,numtopsurface);
        G = embedsurfacesinvolume(G,numsurfaces,numvolume);
        if ~noduplicate
            physicalgroup = 1;
            openboundaryphysicalgroup = 1;
            G = createphysicalcurve(G,[{numembeddedlines},numtoplines],openboundaryphysicalgroup);
            G = createphysicalsurface(G,numsurfaces,physicalgroup);
        end
    end
    if ischarin('recombine',varargin)
        G = recombinesurface(G,{numsurface,numtopsurface});
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
        VIn = min(clH);
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

if ~noduplicate && ...
        (((dim==2 || dim==3 && extrusion) && ~isempty(numembeddedlines)) ...
            || (dim==3 && ~extrusion && ~isempty(numembeddedsurfaces)))
    if (dim==2 && ~isempty(numembeddedlines)) || (dim==3 && ~extrusion && ~isempty(numembeddedsurfaces))
        G = createcrack(G,dim-1,physicalgroup);
    elseif dim==3 && extrusion && ~isempty(numembeddedlines)
        G = createcrack(G,dim-1,physicalgroup,openboundaryphysicalgroup);
    end
    G = remesh(G,dim,varargin{:});
    G = deleteoptfile(G);
    
    [varargout{:}] = gmsh2femobject(indim,G,dim:-1:dim-n+1,varargin{:});
end
