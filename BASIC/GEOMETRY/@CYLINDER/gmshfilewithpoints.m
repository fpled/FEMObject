function varargout = gmshfilewithpoints(C,P,clC,clP,numbercenter,numberpoints,numberembeddedpoints,numbercurves,numbercurveloops,numbersurfaces,numbersurfaceloop,numbervolume,varargin)
% function G = gmshfilewithpoints(C,P,clC,clP,numbercenter,numberpoints,numberembeddedpoints,numbercurves,numbercurveloops,numbersurfaces,numbersurfaceloop,numbervolume)
% C : CYLINDER
% P : POINT
% clC, clP : characteristic length

tol = getfemobjectoptions('tolerancepoint');
angle = C.angle;

if ~iscell(P), P = {P}; end
if nargin<=3 || isempty(clP), clP = clC; end
if isscalar(clP), clP = repmat(clP,1,length(P)); end

% Generate all vertices (base and top)
PC = getvertices(C);
n = numel(PC)/2; % number of points at base/top

if ischarin('extrude',varargin)
    if nargin<=4 || isempty(numbercenter)
        numbercenter = 1; % 1 base center
    end
    if nargin<=5 || isempty(numberpoints)
        numberpoints = 1+(1:n); % n points
    end
    if nargin<=6 || isempty(numberembeddedpoints)
        numberembeddedpoints = 2+5*n+(1:length(P));
    end
    if nargin<=7 || isempty(numbercurves)
        if abs(angle - 2*pi) < tol
            % Full cylinder: closed contour
            numbercurves = 1:n; % n circle arcs
        else
            % Partial cylinder: open arc
            numbercurves = 1:(n+1); % (n-1) circle arcs + 2 radial lines
        end
    end
    if nargin<=8 || isempty(numbercurveloops)
        numbercurveloops = 1; % 1 base circle
    end
else
    if nargin<=4 || isempty(numbercenter)
        numbercenter = 1:2; % 1 base center + 1 top center
    end
    if nargin<=5 || isempty(numberpoints)
        numberpoints = 2+(1:2*n); % n base points + n top points
    end
    if nargin<=6 || isempty(numberembeddedpoints)
        numberembeddedpoints = 2+2*n+(1:length(P));
    end
    if nargin<=7 || isempty(numbercurves)
        if abs(angle - 2*pi) < tol
            % Full cylinder: closed contour
            numbercurves = 1:(3*n); % n base arcs + n top arcs + n verticals
        else
            % Partial cylinder: open arc
            numbercurves = 1:(3*n+3); % (n-1) base arcs + 2 base radials + (n-1) top arcs + 2 top radials + n verticals + 1 center vertical
        end
    end
    if nargin<=8 || isempty(numbercurveloops)
        if abs(angle - 2*pi) < tol
            % Full cylinder: closed contour
            numbercurveloops = 1:(n+2); % 1 base circle + 1 top circle + n lateral faces
        else
            % Partial cylinder: open arc
            numbercurveloops = 1:(n+3); % 1 base circle + 1 top circle + (n-1) lateral faces +  2 radial faces
        end
    end
end
if nargin<=9 || isempty(numbersurfaces), numbersurfaces = numbercurveloops; end
if nargin<=10 || isempty(numbersurfaceloop), numbersurfaceloop = 1; end
if nargin<=11 || isempty(numbervolume), numbervolume = 1; end

G = GMSHFILE();
G = createpoint(G,[C.cx,C.cy,C.cz],clC,numbercenter(1));
if ischarin('extrude',varargin)
    G = createpoints(G,PC(1:n),clC,numberpoints);
    if abs(angle - 2*pi) < tol
        % Full cylinder: closed contour
        G = createcirclecontour(G,numbercenter,numberpoints,numbercurves,numbercurveloops,varargin{:});
    else
        % Partial cylinder: open circle arc with two radials
        G = createcirclearccontour(G,numbercenter,numberpoints,numbercurves,numbercurveloops,varargin{:});
    end
    G = createplanesurface(G,numbercurveloops,numbersurfaces);
    vect = C.h * [C.nx,C.ny,C.nz];
    [G,out] = extrude(G,vect,'Surface',numbersurfaces,varargin{:});
    numbertopsurface = [out,'[0]'];
    numbervolume     = [out,'[1]'];
    nlateralsurfaces = numel(numbercurves);
    numberlateralsurfaces = cell(1,nlateralsurfaces);
    for k=1:nlateralsurfaces
        numberlateralsurfaces{k} = [out, sprintf('[%d]', k+1)];
    end
    numbersurfaces = [{numbersurfaces}, {numbertopsurface}, numberlateralsurfaces{:}];
else
    Ptop = getctop(C);
    G = createpoint(G,Ptop,clC,numbercenter(2));
    G = createpoints(G,PC,clC,numberpoints);
    if abs(angle - 2*pi) < tol
        % Full cylinder: closed contour
        G = createcylindercontour(G,numbercenter,numberpoints,numbercurves,numbercurveloops,numbersurfaces,numbersurfaceloop);
    else
        % Partial cylinder: open circle arc with two radials
        G = createcylinderarccontour(G,numbercenter,numberpoints,numbercurves,numbercurveloops,numbersurfaces,numbersurfaceloop);
    end
    G = createvolume(G,numbersurfaceloop,numbervolume);
end

G = createpoints(G,P,clP,numberembeddedpoints);
G = embedpointsinvolume(G,numberembeddedpoints,numbervolume);

if ischarin('recombine',varargin)
    G = recombinesurface(G);
end

varargout{1} = G;
varargout{2} = numbervolume;
if ischarin('extrude',varargin)
    varargout{3} = numbersurfaces;
end
