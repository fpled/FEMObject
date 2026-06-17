function varargout = gmshThreePointBendingWithSingleEdgeRectangularNotchStructured(L,H,ls,w,a,b,c,t,cl,filename,indim,varargin)
% function varargout = gmshThreePointBendingWithSingleEdgeRectangularNotchStructured(L,H,ls,w,a,b,c,t,cl,filename,indim)
% L : length
% H : height
% ls : location of the supports from the outer edges
% w : flat punch / support width
% a : length/depth of the edge notch
% b : location of the notch centerline from the left outer edge
% c : width of the edge notch
% t : thickness
% cl : characteristic length
% filename : file name (optional)
% indim : space dimension (optional, 2 by default)

recombine = ischarin('recombine',varargin);
varargin = delonlycharin({'recombine'},varargin);

if nargin<11 || isempty(indim)
    indim = 2;
end
if nargin<10 || isempty(filename)
    filename = 'gmsh_three_point_bending_with_single_edge_rectangular_notch_structured';
end
if nargin<8 || isempty(t)
    t = 1;
end
if nargin<6 || isempty(b)
    b = L/2;
end

dim = 2;
tol = getfemobjectoptions('tolerancepoint');
% tol = 1e-12*max([1,L,H]);

xL = b - c/2;
xR = b + c/2;

if ~(L>0 && H>0)
    error('L and H must be strictly positive.');
end
if ~(a>0 && a<H)
    error('The notch length a must satisfy 0 < a < H.');
end
if ~(c>0)
    error('The notch width c must be strictly positive.');
end
if xL <= tol || xR >= L-tol
    error('The notch must satisfy 0 < b-c/2 and b+c/2 < L.');
end
if ~(cl>0)
    error('The characteristic length cl must be strictly positive.');
end

% -------------------------------------------------------------------------
% Global x-grid
% -------------------------------------------------------------------------
% The support and punch points are included as global vertical grid lines.
% This avoids distorted transfinite surfaces.

xSupportLeft  = [];
xSupportRight = [];
xPunch        = [];

if ~isempty(ls) && ~isempty(w) && w>0
    if ls-w/2 < -tol || ls+w/2 > L+tol || ...
       L-ls-w/2 < -tol || L-ls+w/2 > L+tol || ...
       b-w/2 < -tol || b+w/2 > L+tol
        error('Support or punch interval lies outside the domain.');
    end
    xSupportLeft  = [ls-w/2, ls, ls+w/2];
    xSupportRight = [L-ls-w/2, L-ls, L-ls+w/2];
    xPunch        = [b-w/2, b, b+w/2];
end

x = cleanBreakpoints([ ...
    0, L, ...
    xSupportLeft, ...
    xL, b, xR, ...
    xPunch, ...
    xSupportRight],0,L,tol);

% Minimal y-grid: bottom, notch tip level, top.
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
        if isPointInDomain(x(i),y(j))
            numpoints = numpoints + 1;
            P(j,i) = numpoints;
            G = createpoint(G,[x(i),y(j)],cl,P(j,i));
        end
    end
end

% -------------------------------------------------------------------------
% Horizontal curves
% -------------------------------------------------------------------------
Hlines = zeros(ny,nx-1);
numlines = 0;

for j=1:ny
    for i=1:nx-1
        if P(j,i)~=0 && P(j,i+1)~=0 && ...
                isHorizontalSegmentInDomain(x(i),x(i+1),y(j))
            
            numlines = numlines + 1;
            Hlines(j,i) = numlines;
            
            G = createline(G,[P(j,i),P(j,i+1)],Hlines(j,i));
            G = createtransfinitecurve(G,Hlines(j,i),Nx(i)+1);
        end
    end
end

% -------------------------------------------------------------------------
% Vertical curves
% -------------------------------------------------------------------------
Vlines = zeros(ny-1,nx);

for j=1:ny-1
    for i=1:nx
        if P(j,i)~=0 && P(j+1,i)~=0 && ...
                isVerticalSegmentInDomain(x(i),y(j),y(j+1))
            
            numlines = numlines + 1;
            Vlines(j,i) = numlines;
            
            G = createline(G,[P(j,i),P(j+1,i)],Vlines(j,i));
            G = createtransfinitecurve(G,Vlines(j,i),Ny(j)+1);
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
        xc = (x(i)+x(i+1))/2;
        yc = (y(j)+y(j+1))/2;
        
        if isCellInDomain(xc,yc) && ...
                Hlines(j,i)~=0 && Vlines(j,i+1)~=0 && ...
                Hlines(j+1,i)~=0 && Vlines(j,i)~=0
            
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
end

if recombine
    G = recombinesurface(G,numsurfaces);
end

numphysicalsurface = 1;
G = createphysicalsurface(G,numsurfaces,numphysicalsurface);

n = max(nargout,1);
varargout = cell(1,n);
[varargout{:}] = gmsh2femobject(indim,G,dim:-1:dim-n+1,varargin{:});

% -------------------------------------------------------------------------
% Local functions
% -------------------------------------------------------------------------
function rep = isPointInDomain(xp,yp)
    % Remove only the interior of the rectangular notch.
    % Boundary points x=xL, x=xR, y=a are kept.
    rep = ~((xp>xL+tol) && (xp<xR-tol) && (yp<a-tol));
end

function rep = isCellInDomain(xp,yp)
    rep = ~((xp>xL+tol) && (xp<xR-tol) && (yp<a-tol));
end

function rep = isHorizontalSegmentInDomain(x1,x2,yp)
    xm = (x1+x2)/2;
    rep = ~((xm>xL+tol) && (xm<xR-tol) && (yp<a-tol));
end

function rep = isVerticalSegmentInDomain(xp,y1,y2)
    ym = (y1+y2)/2;
    rep = ~((xp>xL+tol) && (xp<xR-tol) && (ym<a-tol));
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