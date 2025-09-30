function varargout = gmshfile_occ(C,cl,numbersurface,varargin)
% function G = gmshfile_occ(C,cl,numbersurface)
% C : CIRCLE
% cl : characteristic length

if nargin<3 || isempty(numbersurface), numbersurface = 1; end

center = [C.cx,C.cy,C.cz];
radius = C.r;

v = [C.vx,C.vy];
n = [C.nx,C.ny,C.nz];
[dir,ang] = calcrotation_direction_angle(C,v,n);

G = GMSHFILE();
G = setfactory(G,"OpenCASCADE");
G = createdisk(G,center,radius,numbersurface);
if abs(ang) > eps
    G = rotate(G,dir,center,ang,'Surface',numbersurface);
end
G = setmeshsize(G,cl);

if ischarin('recombine',varargin)
    G = recombinesurface(G,numbersurface);
end

varargout{1} = G;
varargout{2} = numbersurface;
