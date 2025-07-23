function varargout = gmshfilewithpoints_occ(E,P,clE,clP,numberembeddedpoints,numbersurface,varargin)
% function G = gmshfilewithpoints_occ(E,P,clE,clP,numberembeddedpoints,numbersurface)
% E : ELLIPSE
% P : POINT
% clE, clP : characteristic length

if ~iscell(P), P = {P}; end
if nargin<=3 || isempty(clP), clP = clE; end
if isscalar(clP), clP = repmat(clP,1,length(P)); end

if nargin<=4 || isempty(numberembeddedpoints), numberembeddedpoints = 1+(1:length(P)); end
if nargin<=5 || isempty(numbersurface), numbersurface = 1; end

G = GMSHFILE();
G = setfactory(G,"OpenCASCADE");
G = createdisk(G,[E.cx,E.cy,E.cz],E.a,numbersurface,E.b);
G = setmeshsize(G,clE);

G = createpoints(G,P,clP,numberembeddedpoints);
G = embedpointsinsurface(G,numberembeddedpoints,numbersurface);

if ischarin('recombine',varargin)
    G = recombinesurface(G,numbersurface);
end

varargout{1} = G;
varargout{2} = numbersurface;
