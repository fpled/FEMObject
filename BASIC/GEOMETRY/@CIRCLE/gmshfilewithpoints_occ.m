function varargout = gmshfilewithpoints_occ(C,P,clC,clP,numberembeddedpoints,numbersurface,varargin)
% function G = gmshfilewithpoints_occ(C,P,clC,clP,numberembeddedpoints,numbersurface)
% C : CIRCLE
% P : POINT
% clC, clP : characteristic length

if ~iscell(P), P = {P}; end
if nargin<4 || isempty(clP), clP = clC; end
if isscalar(clP), clP = repmat(clP,1,length(P)); end

if nargin<5 || isempty(numberembeddedpoints), numberembeddedpoints = 1+(1:length(P)); end
if nargin<6 || isempty(numbersurface), numbersurface = 1; end

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
G = setmeshsize(G,clC);

G = createpoints(G,P,clP,numberembeddedpoints);
G = embedpointsinsurface(G,numberembeddedpoints,numbersurface);

if ischarin('recombine',varargin)
    G = recombinesurface(G,numbersurface);
end

varargout{1} = G;
varargout{2} = numbersurface;
