function varargout = gmshfilewithpoints(D,P,clD,clP,numberpoints,numberlines,numberlineloop,numberembeddedpoints,numbersurface,varargin)
% function G = gmshfilewithpoints(D,P,clD,clP,numberpoints,numberlines,numberlineloop,numberembeddedpoints,numbersurface)
% D : QUADRANGLE
% P : POINT
% clD, clP : characteristic length

if ~iscell(P)
    P = {P};
end
if nargin<=3
    clP = clD;
end
if isscalar(clP)
    clP = repmat(clP,1,length(P));
end
if nargin<=4
    numberpoints = 1:4;
    numberlines = 1:4;
    numberlineloop = 5;
    numberembeddedpoints = 4+(1:length(P));
    numbersurface = 1;
elseif nargin==8
    numbersurface = [];
end

G = GMSHFILE();
PD = getvertices(D);
G = createpoints(G,PD(1:4),clD,numberpoints);
G = createcontour(G,numberpoints(1:4),numberlines,numberlineloop);
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
