function varargout = gmshfilewithpoints(C,P,clC,clP,numbercenter,numberpoints,numberlines,numberlineloop,numberembeddedpoints,numbersurface,varargin)
% function G = gmshfilewithpoints(C,P,clC,clP,numbercenter,numberpoints,numberlines,numberlineloop,numberembeddedpoints,numbersurface)
% C : CYLINDER
% P : POINT
% clC, clP : characteristic length

if ~iscell(P)
    P = {P};
end
if nargin<=3
    clP = clC;
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
PC = getvertices(C);
G = createpoint(G,[C.cx,C.cy,C.cz],clC,numbercenter);
G = createpoints(G,PC(1:4),clC,numberpoints);
G = createcirclecontour(G,numbercenter,numberpoints,numberlines,numberlineloop,varargin{:});
if ~isempty(numbersurface)
    G = createplanesurface(G,numberlineloop,numbersurface);
    n = [C.nx,C.ny,C.nz];
    vect = C.h*n;
    G = extrude(G,vect,'Surface',numbersurface,varargin{:});
end
G = createpoints(G,P,clP,numberembeddedpoints);
if ~isempty(numbersurface)
    G = embedpointsinvolume(G,numberembeddedpoints,numbersurface);
    if ischarin('recombine',varargin)
        G = recombinesurface(G);
    end
end

varargout{1} = G;
varargout{2} = numbersurface;
