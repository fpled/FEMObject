function varargout = gmshfilewithpoints_occ(S,P,clS,clP,numberembeddedpoints,numbervolume,varargin)
% function G = gmshfilewithpoints_occ(S,P,clS,clP,numberembeddedpoints,numbervolume)
% S : SPHERE
% P : POINT
% clS, clP : characteristic length

if ~iscell(P), P = {P}; end
if nargin<=3 || isempty(clP), clP = clS; end
if isscalar(clP), clP = repmat(clP,1,length(P)); end

if nargin<=4 || isempty(numberembeddedpoints), numberembeddedpoints = 2+(1:length(P)); end
if nargin<=5 || isempty(numbervolume), numbervolume = 1; end

G = GMSHFILE();
G = setfactory(G,"OpenCASCADE");
G = createsphere(G,[S.cx,S.cy,S.cz],S.r,numbervolume,varargin{:});
G = setmeshsize(G,clS);

G = createpoints(G,P,clP,numberembeddedpoints);
G = embedpointsinvolume(G,numberembeddedpoints,numbervolume);

if ischarin('recombine',varargin)
    G = recombinesurface(G);
end

varargout{1} = G;
varargout{2} = numbervolume;
