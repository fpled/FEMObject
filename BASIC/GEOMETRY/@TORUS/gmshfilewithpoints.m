function varargout = gmshfilewithpoints(T,P,clT,clP,numbercenter,numberpoints,numberembeddedpoints,numbercurves,numbercurveloop,numbersurface,varargin)
% function G = gmshfilewithpoints(C,P,clC,clP,numbercenter,numberpoints,numberembeddedpoints,numbercurves,numbercurveloop,numbersurface)
% C : CYLINDER
% P : POINT
% clC, clP : characteristic length

if ~iscell(P), P = {P}; end
if nargin<=3 || isempty(clP), clP = clT; end
if isscalar(clP) || isempty(clP), clP = repmat(clP,1,length(P)); end

if nargin<=4 || isempty(numbercenter), numbercenter = 1; end
if nargin<=5 || isempty(numberembeddedpoints), numberembeddedpoints = 8+(1:length(P)); end
if nargin<=6 || isempty(numberpoints), numberpoints = 2:5; end
if nargin<=7 || isempty(numbercurves), numbercurves = 1:4; end
if nargin<=8 || isempty(numbercurveloop), numbercurveloop = 1; numbersurface = 1; end
if nargin==9, numbersurface = []; end

G = GMSHFILE();
PC = getvertices(T);
G = createpoint(G,[T.cx,T.cy,T.cz],clT,numbercenter);
G = createpoints(G,PC(1:4),clT,numberpoints);
G = createcirclecontour(G,numbercenter,numberpoints,numbercurves,numbercurveloop,varargin{:});
if ~isempty(numbersurface)
    G = createplanesurface(G,numbercurveloop,numbersurface);
    n = [T.nx,T.ny,T.nz];
    vect = T.h*n;
    G = extrude(G,vect,'Surface',numbersurface,varargin{:});
end

G = createpoints(G,P,clP,numberembeddedpoints);
G = embedpointsinvolume(G,numberembeddedpoints,numbersurface);

if ischarin('recombine',varargin)
    G = recombinesurface(G);
end

varargout{1} = G;
varargout{2} = numbersurface;
