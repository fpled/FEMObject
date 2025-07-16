function varargout = gmshfilewithpoints(E,P,clE,clP,numbercenter,numberpoints,numberlines,numberlineloop,numberembeddedpoints,numbersurface,varargin)
% function G = gmshfilewithpoints(E,P,clE,clP,numbercenter,numberpoints,numberlines,numberlineloop,numberembeddedpoints,numbersurface)
% E : ELLIPSE
% P : POINT
% clE, clP : characteristic length

if ~iscell(P)
    P = {P};
end
if nargin<=3
    clP = clE;
end
if isscalar(clP)
    clP = repmat(clP,1,length(P));
end
if nargin<=4
    numbercenter = 5;
    numberpoints = 1:4;
    numberlines = 1:4;
    numberlineloop = 5;
    numberembeddedpoints = 5+(1:length(P));
    numbersurface = 1;
elseif nargin==9
    numbersurface = [];
end

G = GMSHFILE();
PE = getvertices(E);
G = createpoint(G,[E.cx,E.cy,E.cz],clE,numbercenter);
G = createpoints(G,PE,clE,numberpoints);
G = createellipsecontour(G,numbercenter,numberpoints,numberlines,numberlineloop);
if ~isempty(numbersurface)
    G = createplanesurface(G,numberlineloop,numbersurface);
end
G = createpoints(G,P,clP,numberembeddedpoints);
if ~isempty(numbersurface)
    G = embedpointsinsurface(G,numberembeddedpoints,numbersurface);
    if ischarin('recombine',varargin)
        G = recombinesurface(G,numbersurface);
    end
end

varargout{1} = G;
varargout{2} = numbersurface;
