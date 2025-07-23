function varargout = gmshfilewithpoints_occ(T,P,clT,clP,numberembeddedpoints,numbervolume,varargin)
% function G = gmshfilewithpoints_occ(T,P,clT,clP,numberembeddedpoints,numbervolume)
% T : TORUS
% P = POINT
% clT, clP : characteristic length

if ~iscell(P), P = {P}; end
if nargin<=3 || isempty(clP), clP = clT; end
if isscalar(clP), clP = repmat(clP,1,length(P)); end

if nargin<=4 || isempty(numberembeddedpoints), numberembeddedpoints = 2+(1:length(P)); end
if nargin<=5 || isempty(numbervolume), numbervolume = 1; end

G = GMSHFILE();
G = setfactory(G,"OpenCASCADE");
G = createtorus(G,[T.cx,T.cy,T.cz],[T.r1,T.r2],numbervolume,T.angle);
G = setmeshsize(G,clT);

G = createpoints(G,P,clP,numberembeddedpoints);
G = embedpointsinvolume(G,numberembeddedpoints,numbervolume);

if ischarin('recombine',varargin)
    G = recombinesurface(G);
end

varargout{1} = G;
varargout{2} = numbervolume;
