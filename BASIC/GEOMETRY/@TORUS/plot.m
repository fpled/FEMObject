function varargout = plot(T,varargin)
% function varargout = plot(T,varargin)

tol = getfemobjectoptions('tolerancepoint');

% Major radius, minor radius, and opening angle
r1 = T.r1;
r2 = T.r2;
angle = T.angle;
if isstring(angle), angle = char(angle); end
if ischar(angle),   angle = str2num(lower(angle)); end

isfull = abs(angle - 2*pi) < tol;

% Angular resolution: 200 points for full major circle (angle = 2*pi), scale for partial arc (angle < 2*pi)
%                      50 points for minor circle
npts = getcharin('npts',varargin,max(2,round(200*angle/(2*pi)))); % points along major circle
mpts = getcharin('mpts',varargin,50);                             % points along minor circle

% Major (u) and minor (v) angles
tu = linspace(0,angle,npts+1)'; % parametric major angle
tv = linspace(0,2*pi,mpts+1)';  % parametric minor angle
if isfull, tu(end) = []; end % remove duplicate
tv(end) = []; % remove duplicate
% [U, V] = meshgrid(tu, tv);

% Minor circle locations along the major ring
base_angles = [0, pi/2, pi, 3*pi/2];
minor_angles = base_angles(base_angles < angle + tol); % include only up to angle
% minor_angles = 0; % always include start point at +x
% if angle > pi/2 - tol  , minor_angles(end+1) = pi/2;   end
% if angle > pi - tol    , minor_angles(end+1) = pi;     end
% if angle > 3*pi/2 - tol, minor_angles(end+1) = 3*pi/2; end
if ~isfull && angle > 0
    minor_angles(end+1) = angle; % add endpoint at angle if not a full major circle
end

% Build parametric torus
% X = (r1 + r2 * cos(V)) .* cos(U);
% Y = (r1 + r2 * cos(V)) .* sin(U);
% Z = r2 * sin(V);

% Rotation matrix
v = [T.vx, T.vy];
n = [T.nx, T.ny, T.nz];
R = calcrotation(T,v,n);

% Center
c = [T.cx, T.cy, T.cz];

% Rotate into global frame and translate to center
% nodecoord = [X(:), Y(:), Z(:)] * R + c;

% Reshape back to 2D grids for surf
% sz = size(X);
% X = reshape(nodecoord(:,1), sz);
% Y = reshape(nodecoord(:,2), sz);
% Z = reshape(nodecoord(:,3), sz);

% Connectivity for major circle and major rings
if isfull
    % Full torus
    connec_major = [1:npts,1]; % closed loop
else
    % Partial torus
    connec_major = [(1:npts)' (2:npts+1)']; % polyline
end
% Connectivity for minor circles
connec_minor = [1:mpts,1]; % closed loop

% Plot using patch
options = patchoptions(T.indim,varargin{:});
hs = ishold;
hold on

H = struct();

% Major circle at tube centerline
x0 = r1 * cos(tu);
y0 = r1 * sin(tu);
z0 = zeros(size(tu));
nodecoord_major = [x0, y0, z0] * R + c;
H.major = patch('Faces',connec_major,'Vertices',nodecoord_major,options{:},'LineStyle',':');

% Major rings at the 4 canonical minor angles
H.rings = gobjects(4,1);
for i=1:4
    v = base_angles(i);
    xr = (r1 + r2*cos(v)) * cos(tu);
    yr = (r1 + r2*cos(v)) * sin(tu);
    zr = r2 * sin(v) * ones(size(tu));
    nodecoord_ring = [xr, yr, zr] * R + c;
    H.rings(i) = patch('Faces',connec_major,'Vertices',nodecoord_ring,options{:});
end

% Minor circles and diameters
H.minor = gobjects(numel(minor_angles),1);
H.diameters = gobjects(numel(minor_angles),2);
for i=1:numel(minor_angles)
    % Minor circles
    u = minor_angles(i);
    xc = (r1 + r2*cos(tv)) .* cos(u);
    yc = (r1 + r2*cos(tv)) .* sin(u);
    zc = r2 * sin(tv);
    nodecoord_minor = [xc, yc, zc] * R + c;
    H.minor(i) = patch('Faces',connec_minor,'Vertices',nodecoord_minor,options{:});
    
    % Two orthogonal diameters (tv=0,pi) and (pi/2,3pi/2)
    [~,i1] = min(abs(tv - 0));
    [~,i2] = min(abs(tv - pi));
    [~,i3] = min(abs(tv - pi/2));
    [~,i4] = min(abs(tv - 3*pi/2));
    p1 = nodecoord_minor(i1,:);
    p2 = nodecoord_minor(i2,:);
    p3 = nodecoord_minor(i3,:);
    p4 = nodecoord_minor(i4,:);
    H.diameters(i,1) = patch('Faces',[1,2],'Vertices',[p1; p2],options{:},'LineStyle',':');
    H.diameters(i,2) = patch('Faces',[1,2],'Vertices',[p3; p4],options{:},'LineStyle',':');
end

H = [H.major; H.rings(:); H.minor(:); H.diameters(:)];

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
