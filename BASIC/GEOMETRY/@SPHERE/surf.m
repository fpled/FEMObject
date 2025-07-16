function varargout = surf(S,varargin)
% function varargout = surf(S,varargin)

npts = getcharin('npts',varargin,200); % angular resolution

% Radius
r = S.r;

% Build parametric sphere
[X,Y,Z] = sphere(npts);

% Flatten into vertex list and scale
nodecoord = [X(:), Y(:), Z(:)] * r;

%% Old version
% Rotate around axis n = [nx, ny, nz] by angle of rotation phi =
% atan2(vy, vx) using tangent vector v = [vx, vy] via Rodrigues'
% rotation formula
%% New version
% Twist the XY plane about z = [0, 0, 1] by phi = atan2(vy, vx)
% using tangent vector v = [vx, vy], then tilt from z axis to normal
% vector n = [nx, ny, nz] so that the circle's normal is n regardless
% of v = [vx, vy]
v = [S.vx, S.vy];
n = [S.nx, S.ny, S.nz];
R = calcrotation(S,v,n);

% Translate to center c = [cx, cy, cz]
c = [S.cx, S.cy, S.cz];

% Rotate and translate
nodecoord = nodecoord * R + c;

% Reshape back to 2D grids for surf
sz = size(X);
X = reshape(nodecoord(:,1), sz);
Y = reshape(nodecoord(:,2), sz);
Z = reshape(nodecoord(:,3), sz);

% Plot side using surf
options = patchoptions(S.indim,varargin{:});
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
