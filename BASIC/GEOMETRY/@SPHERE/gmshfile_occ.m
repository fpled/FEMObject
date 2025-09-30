function varargout = gmshfile_occ(S,cl,numbervolume,varargin)
% function G = gmshfile_occ(S,cl,numbervolume)
% S : SPHERE
% cl : characteristic length

if nargin<3 || isempty(numbervolume), numbervolume = 1; end

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
G = setmeshsize(G,cl);

if ischarin('recombine',varargin)
    G = recombinesurface(G);
end

varargout{1} = G;
varargout{2} = numbervolume;
