function varargout = gmshfilewithpoints_occ(E,P,clE,clP,numberembeddedpoints,numbersurface,varargin)
% function G = gmshfilewithpoints_occ(E,P,clE,clP,numberembeddedpoints,numbersurface)
% E : ELLIPSE
% P : POINT
% clE, clP : characteristic length

if ~iscell(P), P = {P}; end
if nargin<4 || isempty(clP), clP = clE; end
if isscalar(clP), clP = repmat(clP,1,length(P)); end

if nargin<5 || isempty(numberembeddedpoints), numberembeddedpoints = 1+(1:length(P)); end
if nargin<6 || isempty(numbersurface), numbersurface = 1; end

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
G = setmeshsize(G,clE);

G = createpoints(G,P,clP,numberembeddedpoints);
G = embedpointsinsurface(G,numberembeddedpoints,numbersurface);

if ischarin('recombine',varargin)
    G = recombinesurface(G,numbersurface);
end

varargout{1} = G;
varargout{2} = numbersurface;
