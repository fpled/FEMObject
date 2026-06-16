function varargout = gmshThreePointBendingWithSingleEdgeCrackStructured(L,H,ls,w,a,b,t,cl,filename,indim,varargin)
% function varargout = gmshThreePointBendingWithSingleEdgeCrackStructured(L,H,ls,w,a,b,t,cl,filename,indim)
% L : length
% H : height
% ls : location of the supports from the outer edges
% w : flat punch / support width
% a : length/depth of the edge crack
% b : location of the edge crack from the left outer edge
% t : thickness
% cl : characteristic length
% filename : file name (optional)
% indim : space dimension (optional, 2 by default)

noduplicate = ischarin('noduplicate',varargin);
recombine   = ischarin('recombine',varargin);

varargin = delonlycharin({'noduplicate','recombine'},varargin);

if nargin<10 || isempty(indim)
    indim = 2;
end
if nargin<9 || isempty(filename)
    filename = 'gmsh_three_point_bending_with_single_edge_crack_structured';
end
if nargin<7 || isempty(t)
    t = 1;
end
if nargin<6 || isempty(b)
    b = L/2;
end

dim = 2;
tol = getfemobjectoptions('tolerancepoint');
% tol = 1e-12*max([1,L,H]);

if ~(L>0 && H>0)
    error('L and H must be strictly positive.');
end
if ~(a>0 && a<H)
    error('The crack length a must satisfy 0 < a < H.');
end
if b <= tol || b >= L-tol
    error('The crack must satisfy 0 < b < L.');
end
if ~(cl>0)
    error('The characteristic length cl must be strictly positive.');
end

% -------------------------------------------------------------------------
% Global x-grid
% -------------------------------------------------------------------------
% The support and punch points are included as global vertical grid lines.
% This preserves the requested boundary points without creating distorted
% transfinite surfaces.

xSupportLeft  = [];
xSupportRight = [];
xPunch        = [];

if ~isempty(ls) && ~isempty(w) && w>0
    xSupportLeft  = [ls-w/2, ls, ls+w/2];
    xSupportRight = [L-ls-w/2, L-ls, L-ls+w/2];
    xPunch        = [b-w/2, b, b+w/2];
end

x = cleanBreakpoints([ ...
    0, L, ...
    xSupportLeft, ...
    b, ...
    xPunch, ...
    xSupportRight],0,L,tol);

% Minimal y-grid: bottom, crack-tip level, top.
y = [0, a, H];

nx = length(x);
ny = length(y);

Nx = max(1,round(diff(x)/cl));
Ny = max(1,round(diff(y)/cl));

% -------------------------------------------------------------------------
% Gmsh file
% -------------------------------------------------------------------------
G = GMSHFILE();

if ischar(filename) && ~isempty(filename)
    G = setfile(G,filename);
end

% -------------------------------------------------------------------------
% Points
% -------------------------------------------------------------------------
P = zeros(ny,nx);
numpoints = 0;

for j=1:ny
    for i=1:nx
        numpoints = numpoints + 1;
        P(j,i) = numpoints;
        G = createpoint(G,[x(i),y(j)],cl,P(j,i));
    end
end

% -------------------------------------------------------------------------
% Horizontal curves
% -------------------------------------------------------------------------
Hlines = zeros(ny,nx-1);
numlines = 0;

for j=1:ny
    for i=1:nx-1
        numlines = numlines + 1;
        Hlines(j,i) = numlines;
        
        G = createline(G,[P(j,i),P(j,i+1)],Hlines(j,i));
        G = createtransfinitecurve(G,Hlines(j,i),Nx(i)+1);
    end
end

% -------------------------------------------------------------------------
% Vertical curves
% -------------------------------------------------------------------------
Vlines = zeros(ny-1,nx);
numlinecrack = [];

for j=1:ny-1
    for i=1:nx
        numlines = numlines + 1;
        Vlines(j,i) = numlines;
        
        G = createline(G,[P(j,i),P(j+1,i)],Vlines(j,i));
        G = createtransfinitecurve(G,Vlines(j,i),Ny(j)+1);
        
        if abs(x(i)-b)<tol && y(j)>=-tol && y(j+1)<=a+tol
            numlinecrack(end+1) = Vlines(j,i);
        end
    end
end

% -------------------------------------------------------------------------
% Surfaces
% -------------------------------------------------------------------------
numlineloop = 0;
numsurface = 0;
numsurfaces = [];

for j=1:ny-1
    for i=1:nx-1
        numlineloop = numlineloop + 1;
        numsurface = numsurface + 1;
        
        numcurves = [Hlines(j,i), ...
                     Vlines(j,i+1), ...
                    -Hlines(j+1,i), ...
                    -Vlines(j,i)];
        
        corners = [P(j,i), ...
                   P(j,i+1), ...
                   P(j+1,i+1), ...
                   P(j+1,i)];
        
        G = createcurveloop(G,numcurves,numlineloop);
        G = createplanesurface(G,numlineloop,numsurface);
        G = createtransfinitesurface(G,numsurface,corners);
        
        numsurfaces(end+1) = numsurface;
    end
end

if recombine
    G = recombinesurface(G,numsurfaces);
end

if ~noduplicate
    openboundaryphysicalgroup = 1;
    physicalgroup = 1;
    
    iCrack = find(abs(x-b)<tol,1,'first');
    jTip   = find(abs(y-a)<tol,1,'first');
    
    if isempty(iCrack) || isempty(jTip) || isempty(numlinecrack)
        error('Internal error: crack line or crack-tip point was not created.');
    end
    
    G = createphysicalpoint(G,P(jTip,iCrack),openboundaryphysicalgroup);
    G = createphysicalcurve(G,numlinecrack,physicalgroup);
end

numphysicalsurface = 1;
G = createphysicalsurface(G,numsurfaces,numphysicalsurface);

n = max(nargout,1);
varargout = cell(1,n);
[varargout{:}] = gmsh2femobject(indim,G,dim:-1:dim-n+1,varargin{:});

if ~noduplicate
    G = createcrack(G,dim-1,physicalgroup,openboundaryphysicalgroup);
    G = remesh(G,dim,varargin{:});
    
    [varargout{:}] = gmsh2femobject(indim,G,dim:-1:dim-n+1,varargin{:});
end

end

function x = cleanBreakpoints(x,xmin,xmax,tol)
x = x(isfinite(x));
x = x(x>=xmin-tol & x<=xmax+tol);
x(x<xmin) = xmin;
x(x>xmax) = xmax;
x(abs(x-xmin)<tol) = xmin;
x(abs(x-xmax)<tol) = xmax;
x = sort(x(:).');

if isempty(x)
    x = [xmin,xmax];
else
    x = x([true,diff(x)>tol]);
end
end
