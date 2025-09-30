function varargout = gmshfilewithpoints(C,P,clC,clP,numbercenter,numberpoints,numberembeddedpoints,numbercurves,numbercurveloops,numbersurfaces,numbersurfaceloop,numbervolume,varargin)
% function G = gmshfilewithpoints(C,P,clC,clP,numbercenter,numberpoints,numberembeddedpoints,numbercurves,numbercurveloops,numbersurfaces,numbersurfaceloop,numbervolume)
% C : CYLINDER
% P : POINT
% clC, clP : characteristic length

extrusion = ischarin('extrude',varargin);
varargin = delonlycharin('extrude',varargin);

tol = getfemobjectoptions('tolerancepoint');
angle = C.angle;

% Normalize and (try to) evaluate angle
if isstring(angle), angle = char(angle); end
if ischar(angle)
    angle_expr = angle;
    angle_val  = str2num(lower(angle));
else
    angle_expr = num2str(angle);
    angle_val  = angle;
end

% Detect full/partial revolution when numeric, or by normalized char if not evaluable
if ~isempty(angle_val)
    isfull = abs(angle_val - 2*pi) < tol;
else
    s = regexprep(lower(angle_expr),'\s+',''); % ignore case and remove spaces (strip whitespace)
    isfull = strcmp(s,'2*pi');
end

if ~iscell(P), P = {P}; end
if nargin<4 || isempty(clP), clP = clC; end
if isscalar(clP), clP = repmat(clP,1,length(P)); end

% Generate all vertices (base and top)
PC = getvertices(C);
n = numel(PC)/2; % number of points at base/top

if extrusion
    if nargin<5 || isempty(numbercenter), numbercenter = 1; end % 1 base center
    if nargin<6 || isempty(numberpoints), numberpoints = numbercenter+(1:n); end % n base points
    % if nargin<7 || isempty(numberembeddedpoints), numberembeddedpoints = max(numberpoints)+(1:length(P)); end
    if nargin<8 || isempty(numbercurves)
        if isfull, numbercurves = 1:n; % n base arcs
        else,      numbercurves = 1:(n+1); % (n-1) base arcs + 2 base radials
        end
    end
    if nargin<9 || isempty(numbercurveloops), numbercurveloops = 1; end % 1 base circle
else
    if nargin<5 || isempty(numbercenter), numbercenter = 1:2; end % 1 base center + 1 top center
    if nargin<6 || isempty(numberpoints), numberpoints = max(numbercenter)+(1:2*n); end % n base points + n top points
    if nargin<7 || isempty(numberembeddedpoints), numberembeddedpoints = max([numbercenter,numberpoints])+(1:length(P)); end
    if nargin<8 || isempty(numbercurves)
        if isfull, numbercurves = 1:(3*n); % n base arcs + n top arcs + n verticals
        else,      numbercurves = 1:(3*n+3); % (n-1) base arcs + 2 base radials + (n-1) top arcs + 2 top radials + n side verticals + 1 center vertical
        end
    end
    if nargin<9 || isempty(numbercurveloops)
        if isfull, numbercurveloops = 1:(n+2); % 1 base circle + 1 top circle + n lateral faces
        else,      numbercurveloops = 1:(n+3); % 1 base circle + 1 top circle + (n-1) lateral faces +  2 radial faces
        end
    end
end
if nargin<10 || isempty(numbersurfaces), numbersurfaces = numbercurveloops; end
if nargin<11 || isempty(numbersurfaceloop), numbersurfaceloop = 1; end
if nargin<12 || isempty(numbervolume), numbervolume = 1; end

center = [C.cx,C.cy,C.cz];

br = @(tag,k) sprintf('%s[%d]', tag, k); % bracket reference helper

G = GMSHFILE();
if extrusion
    % Build base profile as a closed circle contour (full) or open circle arc + two radial lines (partial)
    numbercenter = numbercenter(1);
    numberpoints = numberpoints(1:n);
    if isfull, numbercurves = numbercurves(1:n);
    else,      numbercurves = numbercurves(1:n+1);
    end
    nlatsurfaces = numel(numbercurves);
    numbercurveloop = numbercurveloops(1);
    numbersurface = numbersurfaces(1);
    G = createpoint(G,center,clC,numbercenter);
    G = createpoints(G,PC(1:n),clC,numberpoints);
    if isfull
        % Full cylinder base: closed circle contour
        G = createcirclecontour(G,numbercenter,numberpoints,numbercurves,numbercurveloop);
    else
        % Partial cylinder base: open circle arc + two radial lines
        G = createcirclearccontour(G,numbercenter,numberpoints,numbercurves,numbercurveloop);
    end
    G = createplanesurface(G,numbercurveloop,numbersurface);
    height = C.h;
    n = [C.nx,C.ny,C.nz];
    normal = n / norm(n);
    axis   = height * normal;
    if ischarin('recombine',varargin) && ~ischarin('Layers',varargin)
        numberlayers = max(1, round(height/clC));
        varargin = [varargin, {'Layers',numberlayers}];
    end
    [G,tag] = extrude(G,axis,'Surface',numbersurface,varargin{:});
    numbertopsurface  = br(tag,0); % top surface
    numbervolume      = br(tag,1); % volume
    numberlatsurfaces = arrayfun(@(k) br(tag,k), 1+(1:nlatsurfaces), 'UniformOutput', false); % lateral surfaces
    numbersurfaces    = [{numbersurface, numbertopsurface}, numberlatsurfaces];
    if ischarin('recombine',varargin)
        G = recombinesurface(G,{numbersurface,numbertopsurface});
    end
else
    center_top = getctop(C);
    G = createpoint(G,center,clC,numbercenter(1));
    G = createpoint(G,center_top,clC,numbercenter(2));
    G = createpoints(G,PC,clC,numberpoints);
    if isfull
        % Full cylinder
        G = createcylindercontour(G,numbercenter,numberpoints,numbercurves,numbercurveloops,numbersurfaces,numbersurfaceloop);
    else
        % Partial cylinder
        G = createcylinderarccontour(G,numbercenter,numberpoints,numbercurves,numbercurveloops,numbersurfaces,numbersurfaceloop);
    end
    G = createvolume(G,numbersurfaceloop,numbervolume);
end

if extrusion
    if nargin<7 || isempty(numberembeddedpoints)
        tagembeddedpoint = 'p';
        G = newp(G,tagembeddedpoint);
        numberembeddedpoints = [{tagembeddedpoint}, arrayfun(@(k) sprintf('%s+%d', tagembeddedpoint, k), 1:(length(P)-1), 'UniformOutput', false)];
    end
end
G = createpoints(G,P,clP,numberembeddedpoints);
G = embedpointsinvolume(G,numberembeddedpoints,numbervolume);

if ischarin('recombine',varargin) && ~extrusion
    G = recombinesurface(G);
end

varargout{1} = G;
varargout{2} = numbervolume;
if extrusion
    varargout{3} = numbersurfaces;
end
