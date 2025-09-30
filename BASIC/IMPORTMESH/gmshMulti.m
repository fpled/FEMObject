function varargout = gmshMulti(I,cl,filename,indim,varargin)
% function varargout = gmshMulti(I,cl,filename,indim)
% I : DOMAIN or CIRCLE or ELLIPSE or QUADRANGLE
% cl : characteristic length
% filename : file name (optional)
% indim : space dimension (optional, max(cellfun(@(x) getindim(x),I)) by default)

recombine = ischarin('recombine',varargin);
varargin = delonlycharin('recombine',varargin);

if nargin<4 || isempty(indim)
    indim = max(cellfun(@(x) getindim(x),I));
end

if ~iscell(I)
    I = {I};
end
if isscalar(cl)
    cl = repmat(cl,1,length(I));
end

G = GMSHFILE();
if nargin>=3 && ischar(filename)
    G = setfile(G,filename);
end

numpoints = 0;
numlines = 0;
numlineloop = 0;
numsurface = 0;
if indim==2
    numphysicalsurface = 1;
elseif indim==3
    numvolume = 1;
    numphysicalvolume = 1;
end
for j=1:length(I)
    % numlineloop = [numlineloop,-numlines(1:end-1)];
    if isa(I{j},'DOMAIN') || isa(I{j},'QUADRANGLE')
        if getdim(I{j})==2
            numpoints = numpoints(end)+(1:4);
            numlines = numlines(end)+(1:4);
            numlineloop = numlineloop(end)+1;
            numsurface = numsurface(end)+1;
            GI = gmshfile(I{j},cl(j),numpoints,numlines,numlineloop,numsurface);
        elseif getdim(I{j})==3
            numpoints = numpoints(end)+(1:8);
            numlines = numlines(end)+(1:12);
            numlineloop = numlineloop(end)+(1:6);
            numsurface = numsurface(end)+(1:6);
            numsurfaceloop = numsurfaceloop(end)+1;
            numvolume = numvolume(end)+1;
            GI = gmshfile(I{j},cl(j),numpoints,numlines,numlineloop,numsurface,numsurfaceloop,numvolume);
        end
    elseif isa(I{j},'CIRCLE') || isa(I{j},'ELLIPSE')
        numpoints = numpoints(end)+(1:5);
        numlines = numlines(end)+(1:4);
        numlineloop = numlineloop(end)+1;
        numsurface = numsurface(end)+1;
        GI = gmshfile(I{j},cl(j),numpoints(1),numpoints(2:end),numlines,numlineloop,numsurface);
    elseif isa(I{j},'SPHERE') || isa(I{j},'ELLIPSOID')
        numpoints = numpoints(end)+(1:7);
        numlines = numlines(end)+(1:12);
        numlineloop = numlineloop(end)+(1:8);
        numsurface = numsurface(end)+(1:8);
        numsurfaceloop = numsurfaceloop(end)+1;
        numvolume = numvolume(end)+1;
        GI = gmshfile(I{j},cl(j),numpoints(1),numpoints(2:end),numlines,numlineloop,numsurface,numsurfaceloop,numvolume);
    elseif isa(I{j},'CYLINDER')
        P = getvertices(I{j});
        tol = getfemobjectoptions('tolerancepoint');
        angle = getangle(I{j});
        if isstring(angle), angle = char(angle); end
        if ischar(angle),   angle = str2num(lower(angle)); end
        isfull = abs(angle - 2*pi) < tol;
        n = numel(P)/2; % number of points at base/top
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
        GI = gmshfile(I{j},cl(j),numpoints(1:2),numpoints(3:end),numlines,numlineloop,numsurface,numsurfaceloop,numvolume);
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
        GI = gmshfile(I{j},cl(j),numpoints(1:3+n),numpoints(3*n+1:end),numlines,numlineloop,numsurface,numsurfaceloop,numvolume);
    end
    if indim==2 && getdim(I{j})==2
        numphysicalsurface = numphysicalsurface(end)+1;
        GI = createphysicalsurface(GI,numsurface,numphysicalsurface);
    elseif indim==3 && getdim(I{j})==3
        numphysicalvolume = numphysicalvolume(end)+1;
        GI = createphysicalvolume(GI,numvolume,numphysicalvolume);
    end
    G = G+GI;
end

if recombine
    G = recombinesurface(G);
end

n = max(nargout,1);
varargout = cell(1,n);
dim = max(cellfun(@(x) getdim(x),I));
[varargout{:}] = gmsh2femobject(indim,G,dim:-1:dim-n+1,varargin{:});
