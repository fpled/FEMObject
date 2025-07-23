function varargout = gmshfilewithpoints_occ(C,P,clC,clP,numberembeddedpoints,numbervolume,varargin)
% function G = gmshfilewithpoints_occ(C,P,clC,clP,numberembeddedpoints,numbervolume)
% C : CYLINDER
% P : POINT
% clC, clP : characteristic length

tol = getfemobjectoptions('tolerancepoint');
angle = C.angle;

if ~iscell(P), P = {P}; end
if nargin<=3 || isempty(clP), clP = clC; end
if isscalar(clP), clP = repmat(clP,1,length(P)); end

if nargin<=4 || isempty(numberembeddedpoints), numberembeddedpoints = 6+(1:length(P)); end
if nargin<=5 || isempty(numbervolume), numbervolume = 1; end

G = GMSHFILE();
G = setfactory(G,"OpenCASCADE");
if abs(angle - 2*pi) < tol
    % Full cylinder: closed contour
    G = createcylinder(G,[C.cx,C.cy,C.cz],[C.nx,C.ny,C.nz],C.r,numbervolume);
else
    % Partial cylinder: (npoints-1) circle arcs + 2 radial lines
    G = createcylinder(G,[C.cx,C.cy,C.cz],[C.nx,C.ny,C.nz],C.r,numbervolume,angle);
end
G = setmeshsize(G,clC);

G = createpoints(G,P,clP,numberembeddedpoints);
G = embedpointsinvolume(G,numberembeddedpoints,numbervolume);

if ischarin('recombine',varargin)
    G = recombinesurface(G);
end

varargout{1} = G;
varargout{2} = numbervolume;
