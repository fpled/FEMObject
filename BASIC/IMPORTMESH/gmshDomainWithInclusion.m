function varargout = gmshDomainWithInclusion(D,I,clD,clI,filename,indim,varargin)
% function varargout = gmshDomainWithInclusion(D,I,clD,clI,filename,indim)
% D : DOMAIN or QUADRANGLE
% I : DOMAIN or CIRCLE or ELLIPSE or QUADRANGLE or LINE or POINT
% clD, clI : characteristic lengths
% filename : file name (optional)
% indim : space dimension (optional, getindim(D) by default)

Box       = getcharin('Box',varargin,[]);
extrusion = ischarin('extrude',varargin);

varargin = delcharin('Box',varargin);
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
numembeddedlines = [];
numembeddedsurfaces = [];
if dim==2 || (dim==3 && extrusion)
    numcurves = numlines;
    numsurface = 1;
    numphysicalsurface = 1;
elseif dim==3
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
                GI = gmshfile(I{j},clI(j),numpoints,numlines,numlineloop,numsurface);
                if dim==3
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
            GI = gmshfile(I{j},clI(j),numpoints(1),numpoints(2:end),numlines,numlineloop,numsurface);
            if dim==3
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
                center = getc(I{j});
                GI = GMSHFILE();
                GI = createpoint(GI,center,clI(j),numpoints(1));
                GI = createpoints(GI,P(1:n),clI(j),numpoints(2:end));
                if isfull
                    % Full cylinder base: closed circle contour
                    GI = createcirclecontour(GI,numpoints(1),numpoints(2:end),numlines,numlineloop);
                else
                    % Partial cylinder base: open circle arc + two radial lines
                    GI = createcirclearccontour(GI,numpoints(1),numpoints(2:end),numlines,numlineloop);
                end
                GI = createplanesurface(GI,numlineloop,numsurface);
                numembeddedsurfaces = [numembeddedsurfaces,numsurface];
                I{j} = setdim(I{j},2);
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
                numvolume = numvolume(end)+1;
                GI = gmshfile(I{j},clI(j),numpoints(1:2),numpoints(3:end),numlines,numlineloop,numsurface,numsurfaceloop,numvolume);
            end
        elseif isa(I{j},'TORUS')
            P = getvertices(I{j});
            tol = getfemobjectoptions('tolerancepoint');
            angle = getangle(I{j});
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
            numvolume = numvolume(end)+1;
            GI = gmshfile(I{j},clI(j),numpoints(1:3+n),numpoints(3+n+1:end),numlines,numlineloop,numsurface,numsurfaceloop,numvolume);
        end
        if (dim==2 || (dim==3 && extrusion)) && getdim(I{j})==2
            if dim==2
                if ischarin('recombine',varargin)
                    GI = recombinesurface(GI,numsurface);
                end
                numphysicalsurface = numphysicalsurface(end)+1;
                GI = createphysicalsurface(GI,numsurface,numphysicalsurface);
            end
            numcurves = [numcurves,-numlines];
        elseif dim==3 && getdim(I{j})==3
            numphysicalvolume = numphysicalvolume(end)+1;
            GI = createphysicalvolume(GI,numvolume,numphysicalvolume);
            numsurfaces = [numsurfaces,-numsurface];
        end
    end
    G = G+GI;
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
    if dim==2 || (dim==3 && extrusion)
        G = embedcurvesinsurface(G,numembeddedlines,numsurface);
    elseif dim==3 && ~extrusion
        G = embedcurvesinvolume(G,numembeddedlines,numvolume);
    end
end

if dim==3 && ~extrusion && ~isempty(numembeddedsurfaces)
    G = embedsurfacesinvolume(G,numembeddedsurfaces,numvolume);
end

if dim==3 && extrusion
    % Extrude domain without inclusions
    P = getvertices(D);
    vect = P{5}-P{1};
    if ischarin('recombine',varargin) && ~ischarin('Layers',varargin)
        numlayers = max(1, round(vect(3)/mean([clD,clI(:)'])));
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
            tag = ['outline' num2str(i)];
            G = extrude(G,vect,'Curve',numembeddedlines(i),tag,varargin{:});
            numtoplines{i} = br(tag,0); % top line
            numsurfaces{i} = br(tag,1); % surface
            % numtoppoints{i} = ['p' num2str(i)];
            % G = pointsofcurve(G,numtoplines{i},numtoppoints{i});
        end
        G = embedcurvesinsurface(G,numtoplines,numtopsurface);
        G = embedsurfacesinvolume(G,numsurfaces,numvolume);
    end
    % Extrude embedded surfaces
    if ~isempty(numembeddedsurfaces)
        N = numel(numembeddedsurfaces);
        numtopsurfaces = cell(1,N);
        numvolumes     = cell(1,N);
        for i=1:N
            tag = ['outsurf' num2str(i)];
            G = extrude(G,vect,'Surface',numembeddedsurfaces(i),tag,varargin{:});
            numtopsurfaces{i} = br(tag,0); % top surface
            numvolumes{i}     = br(tag,1); % volume
        end
        numvolume = [{numvolume},numvolumes];
        if ischarin('recombine',varargin)
            G = recombinesurface(G,[{numembeddedsurfaces},numtopsurfaces]);
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
        VIn = min(clI);
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
