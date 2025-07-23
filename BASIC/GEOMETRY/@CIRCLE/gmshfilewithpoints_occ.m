function varargout = gmshfilewithpoints_occ(C,P,clC,clP,numberembeddedpoints,numbersurface,varargin)
% function G = gmshfilewithpoints_occ(C,P,clC,clP,numberembeddedpoints,numbersurface)
% C : CIRCLE
% P : POINT
% clC, clP : characteristic length

if ~iscell(P), P = {P}; end
if nargin<=3 || isempty(clP), clP = clC; end
if isscalar(clP), clP = repmat(clP,1,length(P)); end

if nargin<=4 || isempty(numberembeddedpoints), numberembeddedpoints = 1+(1:length(P)); end
if nargin<=5 || isempty(numbersurface), numbersurface = 1; end

G = GMSHFILE();
G = setfactory(G,"OpenCASCADE");
G = createdisk(G,[C.cx,C.cy,C.cz],C.r,numbersurface);
G = setmeshsize(G,clC);

G = createpoints(G,P,clP,numberembeddedpoints);
G = embedpointsinsurface(G,numberembeddedpoints,numbersurface);

if ischarin('recombine',varargin)
    G = recombinesurface(G,numbersurface);
end

varargout{1} = G;
varargout{2} = numbersurface;
