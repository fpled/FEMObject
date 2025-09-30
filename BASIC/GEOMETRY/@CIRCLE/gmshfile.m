function varargout = gmshfile(C,cl,numbercenter,numberpoints,numbercurves,numbercurveloop,numbersurface,varargin)
% function G = gmshfile(C,cl,numbercenter,numberpoints,numbercurves,numbercurveloop,numbersurface)
% C : CIRCLE
% cl : characteristic length

if nargin<3 || isempty(numbercenter), numbercenter = 1; end
if nargin<4 || isempty(numberpoints), numberpoints = numbercenter+(1:4); end
if nargin<5 || isempty(numbercurves), numbercurves = 1:4; end
if nargin<6 || isempty(numbercurveloop), numbercurveloop = 1; end
if nargin<7, numbersurface = []; elseif isempty(numbersurface), numbersurface = 1; end

center = [C.cx,C.cy,C.cz];
P = getvertices(C);

G = GMSHFILE();
G = createpoint(G,center,cl,numbercenter);
G = createpoints(G,P,cl,numberpoints);
G = createcirclecontour(G,numbercenter,numberpoints,numbercurves,numbercurveloop);
if ~isempty(numbersurface)
    G = createplanesurface(G,numbercurveloop,numbersurface);
    if ischarin('recombine',varargin)
        G = recombinesurface(G,numbersurface);
    end
end

varargout{1} = G;
varargout{2} = numbersurface;
