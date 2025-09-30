function varargout = gmshfile(D,cl,numberpoints,numberlines,numberlineloop,numbersurface,varargin)
% function G = gmshfile(D,cl,numberpoints,numberlines,numberlineloop,numbersurface)
% D : QUADRANGLE
% cl : characteristic length

if nargin<3 || isempty(numberpoints), numberpoints = 1:4; end
if nargin<4 || isempty(numberlines), numberlines = 1:4; end
if nargin<5 || isempty(numberlineloop), numberlineloop = 1; end
if nargin<6, numbersurface = []; elseif isempty(numbersurface), numbersurface = 1; end

P = getvertices(D);

G = GMSHFILE();
G = createpoints(G,P,cl,numberpoints);
G = createcontour(G,numberpoints,numberlines,numberlineloop);
if ~isempty(numbersurface)
    G = createplanesurface(G,numberlineloop,numbersurface);
    if ischarin('recombine',varargin)
        G = recombinesurface(G,numbersurface);
    end
end

varargout{1} = G;
varargout{2} = numbersurface;
