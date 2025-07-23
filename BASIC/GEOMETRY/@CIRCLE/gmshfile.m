function varargout = gmshfile(C,cl,numbercenter,numberpoints,numbercurves,numbercurveloop,numbersurface,varargin)
% function G = gmshfile(C,cl,numbercenter,numberpoints,numbercurves,numbercurveloop,numbersurface)
% C : CIRCLE
% cl : characteristic length

if nargin<=2 || isempty(numbercenter), numbercenter = 1; end
if nargin<=3 || isempty(numberpoints), numberpoints = 2:5; end
if nargin<=4 || isempty(numbercurves), numbercurves = 1:4; end
if nargin<=5 || isempty(numbercurveloop), numbercurveloop = 1; numbersurface = 1; end
if nargin==6, numbersurface = []; end

G = GMSHFILE();
P = getvertices(C);
G = createpoint(G,[C.cx,C.cy,C.cz],cl,numbercenter);
G = createpoints(G,P,cl,numberpoints);
G = createcirclecontour(G,numbercenter,numberpoints,numbercurves,numbercurveloop,varargin{:});
if ~isempty(numbersurface)
    G = createplanesurface(G,numbercurveloop,numbersurface);
    if ischarin('recombine',varargin)
        G = recombinesurface(G,numbersurface);
    end
end

varargout{1} = G;
varargout{2} = numbersurface;
