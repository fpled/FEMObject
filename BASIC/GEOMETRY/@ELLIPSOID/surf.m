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

%% Old version
% Rotate around axis n = [nx, ny, nz] by angle of rotation phi =
% atan2(vy, vx) using tangent vector v = [vx, vy] via Rodrigues'
% rotation formula
%% New version
% Twist the XY plane about z = [0, 0, 1] by phi = atan2(vy, vx)
% using tangent vector v = [vx, vy], then tilt from z axis to normal
% vector n = [nx, ny, nz] so that the circle's normal is n regardless
% of v = [vx, vy]
v = [E.vx, E.vy];
n = [E.nx, E.ny, E.nz];
R = calcrotation(E,v,n);

% Translate to center ct = [cx, cy, cz]
ct = [E.cx, E.cy, E.cz];

% Rotate and translate
nodecoord = nodecoord * R + ct;

% Reshape back to 2D grids for surf
sz = size(X);
X = reshape(nodecoord(:,1), sz);
Y = reshape(nodecoord(:,2), sz);
Z = reshape(nodecoord(:,3), sz);

% Plot side using surf
options = patchoptions(E.indim,varargin{:});
H = surf(X, Y, Z, options{:});

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

if nargout>=1
    varargout{1} = H;
end
