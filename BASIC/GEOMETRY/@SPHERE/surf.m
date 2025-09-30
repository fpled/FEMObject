function varargout = surf(S,varargin)
% function varargout = surf(S,varargin)

npts = getcharin('npts',varargin,200); % angular resolution

% Radius
r = S.r;

% Build parametric sphere
[X,Y,Z] = sphere(npts);

% Flatten into vertex list and scale
nodecoord = [X(:), Y(:), Z(:)] * r;

% Rotation matrix
v = [S.vx, S.vy];
n = [S.nx, S.ny, S.nz];
R = calcrotation(S,v,n);

% Center
c = [S.cx, S.cy, S.cz];

% Rotate into global frame and translate to center
nodecoord = nodecoord * R + c;

% Reshape back to 2D grids for surf
sz = size(X);
X = reshape(nodecoord(:,1), sz);
Y = reshape(nodecoord(:,2), sz);
Z = reshape(nodecoord(:,3), sz);

% Plot using surf
options = patchoptions(S.indim,varargin{:});
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
