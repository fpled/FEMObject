function varargout = plot(C,varargin)
% function varargout = plot(C,varargin)

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
t = linspace(0,angle,npts+1)'; % parametric angle
if isfull, t(end) = []; end % remove duplicate

% Build parametric cylinder
if isfull
    % Full cylinder
    [X,Y,Z] = cylinder(r,npts);
    Z = Z * h;
    nodecoord_base = [X(1,:)', Y(1,:)', Z(1,:)'];
    nodecoord_top  = [X(2,:)', Y(2,:)', Z(2,:)'];
else
    % Partial cylinder
    x = r * cos(t);
    y = r * sin(t);
    z_bot = zeros(size(t));
    z_top = h * ones(size(t));
    nodecoord_base = [x, y, z_bot];
    nodecoord_top  = [x, y, z_top];
end

% Rotation matrix
v = [C.vx, C.vy];
n = [C.nx, C.ny, C.nz];
R = calcrotation(C,v,n);

% Center
c = [C.cx, C.cy, C.cz];

% Rotate into global frame and translate to center
nodecoord_base = nodecoord_base * R + c;
nodecoord_top  = nodecoord_top * R + c;
center_base = [0, 0, 0] * R + c;
center_top  = [0, 0, h] * R + c;

% Plot using patch
options = patchoptions(C.indim,varargin{:});
hs = ishold;
hold on

H = struct();

% Base and top circles
H.circles = gobjects(2,1);
if isfull
    % Full circle: use patch for closed loop
    connec = [1:npts,1]; % connectivity for a closed loop
    H.circles(1) = patch('Faces',connec,'Vertices',nodecoord_base,options{:}); % base circle
    H.circles(2) = patch('Faces',connec,'Vertices',nodecoord_top,options{:}); % top circle
else
    % Partial cylinder: create a sector (fan) for open circle arc by including center point
    nodecoord_base_fan = [center_base; nodecoord_base]; % add base center
    nodecoord_top_fan  = [center_top;  nodecoord_top];  % add top center
    % connec_fan = [1:(length(t)+1),1]; % connectivity for a closed loop
    connec_fan = [1:(npts+2),1]; % connectivity for a closed loop
    H.circles(1) = patch('Faces',connec_fan,'Vertices',nodecoord_base_fan,options{:});
    H.circles(2) = patch('Faces',connec_fan,'Vertices',nodecoord_top_fan,options{:});
end

% Vertical spoke along the central axis (center base to center top)
H.axis = patch('Faces',[1 2],'Vertices',[center_base; center_top],options{:},'LineStyle','-.');

% Vertical spoke lines and radial lines
spoke_angles = [0, pi/2, pi, 3*pi/2];
spoke_idx = 1; % start
for k=2:4
    if angle >= spoke_angles(k)-tol
        [~,i] = min(abs(t - spoke_angles(k)));
        spoke_idx = [spoke_idx, i];
    end
end
if ~isfull
    spoke_idx = unique([spoke_idx, length(t)]); % also include endpoint for open circle arc
end
spoke_idx = unique(spoke_idx);

H.spokes  = gobjects(length(spoke_idx),1);
H.radials = gobjects(length(spoke_idx),2);
for k=1:length(spoke_idx)
    i = spoke_idx(k);
    H.spokes(k)    = patch('Faces',[1 2],'Vertices',[nodecoord_base(i,:); nodecoord_top(i,:)],options{:});
    H.radials(k,1) = patch('Faces',[1 2],'Vertices',[center_base; nodecoord_base(i,:)],options{:},'LineStyle',':');
    H.radials(k,2) = patch('Faces',[1 2],'Vertices',[center_top;  nodecoord_top(i,:)],options{:},'LineStyle',':');
end

H = [H.circles(:); H.axis; H.spokes(:); H.radials(:)];

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
