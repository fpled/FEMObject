function varargout = gmshfilewithpoints(D,P,clD,clP,numberpoints,numberembeddedpoints,numberlines,numberlineloop,numbersurface,varargin)
% function G = gmshfilewithpoints(D,P,clD,clP,numberpoints,numberembeddedpoints,numberlines,numberlineloop,numbersurface)
% D : QUADRANGLE
% P : POINT
% clD, clP : characteristic length

if ~iscell(P), P = {P}; end
if nargin<4 || isempty(clP), clP = clD; end
if isscalar(clP), clP = repmat(clP,1,length(P)); end

if nargin<5 || isempty(numberpoints), numberpoints = 1:4; end
if nargin<6 || isempty(numberembeddedpoints), numberembeddedpoints = max(numberpoints)+(1:length(P)); end
if nargin<7 || isempty(numberlines), numberlines = 1:4; end
if nargin<8 || isempty(numberlineloop), numberlineloop = 1; end
if nargin<9 || isempty(numbersurface), numbersurface = 1; end

PD = getvertices(D);

G = GMSHFILE();
G = createpoints(G,PD,clD,numberpoints);
G = createcontour(G,numberpoints,numberlines,numberlineloop);
G = createplanesurface(G,numberlineloop,numbersurface);

G = createpoints(G,P,clP,numberembeddedpoints);
G = embedpointsinsurface(G,numberembeddedpoints,numbersurface);

if ischarin('recombine',varargin)
    G = recombinesurface(G,numbersurface);
end

varargout{1} = G;
varargout{2} = numbersurface;
