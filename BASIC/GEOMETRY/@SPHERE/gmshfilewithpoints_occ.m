function varargout = gmshfilewithpoints_occ(S,P,clS,clP,numberembeddedpoints,numbervolume,varargin)
% function G = gmshfilewithpoints_occ(S,P,clS,clP,numberembeddedpoints,numbervolume)
% S : SPHERE
% P : POINT
% clS, clP : characteristic length

if ~iscell(P), P = {P}; end
if nargin<4 || isempty(clP), clP = clS; end
if isscalar(clP), clP = repmat(clP,1,length(P)); end

if nargin<5 || isempty(numberembeddedpoints), numberembeddedpoints = 2+(1:length(P)); end
if nargin<6 || isempty(numbervolume), numbervolume = 1; end

center = [S.cx,S.cy,S.cz];
radius = S.r;

v = [S.vx,S.vy];
n = [S.nx,S.ny,S.nz];
[dir,ang] = calcrotation_direction_angle(S,v,n);

G = GMSHFILE();
G = setfactory(G,"OpenCASCADE");
G = createsphere(G,center,radius,numbervolume);
if abs(ang) > eps
    G = rotate(G,dir,center,ang,'Volume',numbervolume);
end
G = setmeshsize(G,clS);

G = createpoints(G,P,clP,numberembeddedpoints);
G = embedpointsinvolume(G,numberembeddedpoints,numbervolume);

if ischarin('recombine',varargin)
    G = recombinesurface(G);
end

varargout{1} = G;
varargout{2} = numbervolume;
