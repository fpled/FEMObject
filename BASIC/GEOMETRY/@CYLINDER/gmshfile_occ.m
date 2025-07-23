function varargout = gmshfile_occ(C,cl,numbervolume,varargin)
% function G = gmshfile_occ(C,cl,numbervolume)
% C : CYLINDER
% cl : characteristic length

tol = getfemobjectoptions('tolerancepoint');
angle = C.angle;

if nargin<=2 || isempty(numbervolume), numbervolume = 1; end

G = GMSHFILE();
G = setfactory(G,"OpenCASCADE");
if abs(angle - 2*pi) < tol
    % Full cylinder: closed contour
    G = createcylinder(G,[C.cx,C.cy,C.cz],[C.nx,C.ny,C.nz],C.r,numbervolume);
else
    % Partial cylinder: (n-1) circle arcs + 2 radial lines
    G = createcylinder(G,[C.cx,C.cy,C.cz],[C.nx,C.ny,C.nz],C.r,numbervolume,angle);
end
G = setmeshsize(G,cl);

if ischarin('recombine',varargin)
    G = recombinesurface(G);
end

varargout{1} = G;
varargout{2} = numbervolume;
