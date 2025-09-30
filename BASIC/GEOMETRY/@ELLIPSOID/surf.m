function varargout = surf(E,varargin)
% function varargout = surf(E,varargin)

npts = getcharin('npts',varargin,200); % angular resolution

% Semi-axes
a = E.a;
b = E.b;
c = E.c;

% Build parametric ellipsoid
[X,Y,Z] = ellipsoid(E.cx,E.cy,E.cz,a,b,c,npts);

% Flatten into vertex list
nodecoord = [X(:), Y(:), Z(:)];

% Rotation matrix
v = [E.vx, E.vy];
n = [E.nx, E.ny, E.nz];
R = calcrotation(E,v,n);

% Center
ct = [E.cx, E.cy, E.cz];

% Rotate into global frame and translate to center
nodecoord = nodecoord * R + ct;

% Reshape back to 2D grids for surf
sz = size(X);
X = reshape(nodecoord(:,1), sz);
Y = reshape(nodecoord(:,2), sz);
Z = reshape(nodecoord(:,3), sz);

% Plot side using surf
options = patchoptions(E.indim,varargin{:});
hs = ishold;
hold on

H = surf(X, Y, Z, options{:});

if ~hs
    hold off
end

axis image

% Optional view or camera controls
numview = getcharin('view',varargin);
up_vector = getcharin('camup',varargin);
camera_position = getcharin('campos',varargin);
if ~isempty(numview)
    view(numview)
else
    view(3)
end
if ~isempty(up_vector)
    camup(up_vector)
end
if ~isempty(camera_position)
    campos(camera_position);
end

if nargout
    varargout{1} = H;
end
