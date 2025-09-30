function varargout = surf(T,varargin)
% function varargout = surf(T,varargin)

% Major radius, minor radius, and opening angle
r1 = T.r1;
r2 = T.r2;
angle = T.angle;
if isstring(angle), angle = char(angle); end
if ischar(angle),   angle = str2num(lower(angle)); end

% Angular resolution: 200 points for full major circle (angle = 2*pi), scale for partial arc (angle < 2*pi)
%                     100 points for minor circle
npts = getcharin('npts',varargin,max(2,round(200*angle/(2*pi)))); % points along major circle
mpts = getcharin('mpts',varargin,50);                            % points along minor circle

% Major (u) and minor (v) angles
tu = linspace(0,angle,npts+1)'; % parametric major angle
tv = linspace(0,2*pi,mpts+1)';  % parametric minor angle
[U, V] = meshgrid(tu, tv);

% Build parametric torus
X = (r1 + r2 * cos(V)) .* cos(U);
Y = (r1 + r2 * cos(V)) .* sin(U);
Z = r2 * sin(V);

% Rotation matrix
v = [T.vx, T.vy];
n = [T.nx, T.ny, T.nz];
R = calcrotation(T,v,n);

% Center
c = [T.cx, T.cy, T.cz];

% Rotate into global frame and translate to center
nodecoord = [X(:), Y(:), Z(:)] * R + c;

% Reshape back to 2D grids for surf
sz = size(X);
X = reshape(nodecoord(:,1), sz);
Y = reshape(nodecoord(:,2), sz);
Z = reshape(nodecoord(:,3), sz);

% Plot using surf
options = patchoptions(T.indim,varargin{:});
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
