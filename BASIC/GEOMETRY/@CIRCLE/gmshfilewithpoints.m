function varargout = gmshfilewithpoints(C,P,clC,clP,numbercenter,numberpoints,numberembeddedpoints,numbercurves,numbercurveloop,numbersurface,varargin)
% function G = gmshfilewithpoints(C,P,clC,clP,numbercenter,numberpoints,numberembeddedpoints,numbercurves,numbercurveloop,numbersurface)
% C : CIRCLE
% P : POINT
% clC, clP : characteristic length

if ~iscell(P), P = {P}; end
if nargin<=3 || isempty(clP), clP = clC; end
if isscalar(clP), clP = repmat(clP,1,length(P)); end

if nargin<=4 || isempty(numbercenter), numbercenter = 1; end
if nargin<=5 || isempty(numberpoints), numberpoints = 2:5; end
if nargin<=6 || isempty(numberembeddedpoints), numberembeddedpoints = numberpoints(end)+(1:length(P)); end
if nargin<=7 || isempty(numbercurves), numbercurves = 1:4; end
if nargin<=8 || isempty(numbercurveloop), numbercurveloop = 1; end
if nargin<=9 || isempty(numbersurface), numbersurface = 1; end

G = GMSHFILE();
PC = getvertices(C);
G = createpoint(G,[C.cx,C.cy,C.cz],clC,numbercenter);
G = createpoints(G,PC,clC,numberpoints);
G = createcirclecontour(G,numbercenter,numberpoints,numbercurves,numbercurveloop,varargin{:});
G = createplanesurface(G,numbercurveloop,numbersurface);

G = createpoints(G,P,clP,numberembeddedpoints);
G = embedpointsinsurface(G,numberembeddedpoints,numbersurface);

if ischarin('recombine',varargin)
    G = recombinesurface(G,numbersurface);
end

varargout{1} = G;
varargout{2} = numbersurface;
