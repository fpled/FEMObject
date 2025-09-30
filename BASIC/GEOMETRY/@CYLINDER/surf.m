function varargout = surf(C,varargin)
% function varargout = surf(C,varargin)

tol = getfemobjectoptions('tolerancepoint');

% Radius, height, and opening angle
r = C.r;
h = C.h;
angle = C.angle;
if isstring(angle), angle = char(angle); end
if ischar(angle),   angle = str2num(lower(angle)); end

isfull = abs(angle - 2*pi) < tol;

% Angular resolution: 200 points for full cylinder (angle = 2*pi), scale for partial cylinder (angle < 2*pi)
npts = getcharin('npts',varargin,max(2,round(200*angle/(2*pi))));

% Build parametric cylinder
if isfull
    % Full cylinder
    [X,Y,Z] = cylinder(r,npts);
    Z = Z * h;
else
    % Partial cylinder
    t = linspace(0,angle,npts+1)'; % parametric angle
    z = [0, h];
    [T,Z] = meshgrid(t,z);
    X = r * cos(T);
    Y = r * sin(T);
end

% Flatten into vertex list
nodecoord = [X(:), Y(:), Z(:)];

% Rotation matrix
v = [C.vx, C.vy];
n = [C.nx, C.ny, C.nz];
R = calcrotation(C,v,n);

% Center
c = [C.cx, C.cy, C.cz];

% Rotate into global frame and translate to center
nodecoord = nodecoord * R + c;

% Reshape back to 2D grids for surf
sz = size(X);
X = reshape(nodecoord(:,1), sz);
Y = reshape(nodecoord(:,2), sz);
Z = reshape(nodecoord(:,3), sz);

% Plot using surf
options = patchoptions(C.indim,varargin{:});
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
