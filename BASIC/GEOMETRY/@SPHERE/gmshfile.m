function varargout = gmshfile(S,cl,numbercenter,numberpoints,numberlines,numberlineloop,numbersurface,varargin)
% function G = gmshfile(S,cl,numbercenter,numberpoints,numberlines,numberlineloop,numbersurface)
% S : SPHERE
% cl : characteristic length

if nargin<=2
    numbercenter = 5;
    numberpoints = 1:4;
    numberlines = 1:4;
    numberlineloop = 5;
    numbersurface = 1;
elseif nargin==6
    numbersurface = [];
end

G = GMSHFILE();
P = getvertices(S);
G = createpoint(G,[S.cx,S.cy,S.cz],cl,numbercenter);
G = createpoints(G,P,cl,numberpoints);
G = createcirclecontour(G,numbercenter,numberpoints,numberlines,numberlineloop,varargin{:});
if ~isempty(numbersurface)
    G = createplanesurface(G,numberlineloop,numbersurface);
    if ischarin('recombine',varargin)
        G = recombinesurface(G,numbersurface);
    end
end

varargout{1} = G;
varargout{2} = numbersurface;
