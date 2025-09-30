function varargout = gmshfile_occ(E,cl,numbersurface,varargin)
% function G = gmshfile_occ(E,cl,numbersurface)
% E : ELLIPSE
% cl : characteristic length

if nargin<3 || isempty(numbersurface), numbersurface = 1; end

center  = [E.cx,E.cy,E.cz];
radius  = E.a;
radiusy = E.b;

v = [E.vx,E.vy];
n = [E.nx,E.ny,E.nz];
[dir,ang] = calcrotation_direction_angle(E,v,n);

G = GMSHFILE();
G = setfactory(G,"OpenCASCADE");
G = createdisk(G,center,radius,numbersurface,radiusy);
if abs(ang) > eps
    G = rotate(G,dir,center,ang,'Surface',numbersurface);
end
G = setmeshsize(G,cl);

if ischarin('recombine',varargin)
    G = recombinesurface(G,numbersurface);
end

varargout{1} = G;
varargout{2} = numbersurface;
