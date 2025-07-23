function varargout = gmshfile(C,cl,numbercenter,numberpoints,numbercurves,numbercurveloops,numbersurfaces,numbersurfaceloop,numbervolume,varargin)
% function G = gmshfile(C,cl,numbercenter,numberpoints,numbercurves,numbercurveloops,numbersurfaces,numbersurfaceloop,numbervolume)
% C : CYLINDER
% cl : characteristic length

tol = getfemobjectoptions('tolerancepoint');
angle = C.angle;

% Generate all vertices (base and top)
P = getvertices(C);
n = numel(P)/2; % number of points at base/top

if ischarin('extrude',varargin)
    if nargin<=2 || isempty(numbercenter)
        numbercenter = 1; % 1 base center
    end
    if nargin<=3 || isempty(numberpoints)
        numberpoints = 1+(1:n); % n points
    end
    if nargin<=4 || isempty(numbercurves)
        if abs(angle - 2*pi) < tol
            % Full cylinder: closed contour
            numbercurves = 1:n; % n circle arcs
        else
            % Partial cylinder: open arc
            numbercurves = 1:(n+1); % (n-1) circle arcs + 2 radial lines
        end
    end
    if nargin<=5 || isempty(numbercurveloops)
        numbercurveloops = 1; % 1 base circle
    end
else
    if nargin<=2 || isempty(numbercenter)
        numbercenter = 1:2; % 1 base center + 1 top center
    end
    if nargin<=3 || isempty(numberpoints)
        numberpoints = 2+(1:2*n); % n base points + n top points
    end
    if nargin<=4 || isempty(numbercurves)
        if abs(angle - 2*pi) < tol
            % Full cylinder: closed contour
            numbercurves = 1:(3*n); % n base arcs + n top arcs + n verticals
        else
            % Partial cylinder: open arc
            numbercurves = 1:(3*n+3); % (n-1) base arcs + 2 base radials + (n-1) top arcs + 2 top radials + n verticals + 1 center vertical
        end
    end
    if nargin<=5 || isempty(numbercurveloops)
        if abs(angle - 2*pi) < tol
            % Full cylinder: closed contour
            numbercurveloops = 1:(n+2); % 1 base circle + 1 top circle + n lateral faces
        else
            % Partial cylinder: open arc
            numbercurveloops = 1:(n+3); % 1 base circle + 1 top circle + (n-1) lateral faces +  2 radial faces
        end
    end
end
if nargin<=6 || isempty(numbersurfaces), numbersurfaces = numbercurveloops; end
if nargin<=7, numbersurfaceloop = 1; numbervolume = 1; end
if nargin==8, numbervolume = []; end

G = GMSHFILE();
G = createpoint(G,[C.cx,C.cy,C.cz],cl,numbercenter(1));
if ischarin('extrude',varargin)
    G = createpoints(G,P(1:n),cl,numberpoints);
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
    G = createpoint(G,Ptop,cl,numbercenter(2));
    G = createpoints(G,P,cl,numberpoints);
    if abs(angle - 2*pi) < tol
        % Full cylinder: closed contour
        G = createcylindercontour(G,numbercenter,numberpoints,numbercurves,numbercurveloops,numbersurfaces,numbersurfaceloop);
    else
        % Partial cylinder: open circle arc with two radials
        G = createcylinderarccontour(G,numbercenter,numberpoints,numbercurves,numbercurveloops,numbersurfaces,numbersurfaceloop);
    end
    if ~isempty(numbervolume)
        G = createvolume(G,numbersurfaceloop,numbervolume);
    end
end

if ischarin('recombine',varargin)
    G = recombinesurface(G);
end

varargout{1} = G;
varargout{2} = numbervolume;
if ischarin('extrude',varargin)
    varargout{3} = numbersurfaces;
end
