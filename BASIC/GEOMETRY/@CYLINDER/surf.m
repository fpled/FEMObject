function varargout = surf(C,varargin)
% function varargout = surf(C,varargin)

% Radius, height and opening angle
r = C.r;
h = C.h;
angle = C.angle;

% Angular resolution: 200 points for full circle (angle = 2*pi), scale for partial arc
npts = getcharin('npts',varargin,max(2,round(200*angle/(2*pi))));

% Build parametric cylinder
if angle==2*pi
    [X,Y,Z] = cylinder(r,npts);
    Z = Z * h;
else
    t = linspace(0,angle,npts+1)';  % parametric angle
    z = [0, h];
    [T,Z] = meshgrid(t,z);
    X = r * cos(T);
    Y = r * sin(T);
end

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
v = [C.vx, C.vy];
n = [C.nx, C.ny, C.nz];
R = calcrotation(C,v,n);

% Translate to center c = [cx, cy, cz]
c = [C.cx, C.cy, C.cz];

% Rotate and translate
nodecoord = nodecoord * R + c;

% Reshape back to 2D grids for surf
sz = size(X);
X = reshape(nodecoord(:,1), sz);
Y = reshape(nodecoord(:,2), sz);
Z = reshape(nodecoord(:,3), sz);

% Plot side using surf
options = patchoptions(C.indim,varargin{:});
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
