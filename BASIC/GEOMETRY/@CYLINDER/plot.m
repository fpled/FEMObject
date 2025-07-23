function varargout = plot(C,varargin)
% function varargout = plot(C,varargin)

% Radius, height and opening angle
r = C.r;
h = C.h;
angle = C.angle;

% Angular resolution: 200 points for full circle (angle = 2*pi), scale for partial arc
npts = getcharin('npts',varargin,max(2,round(200*angle/(2*pi))));
t = linspace(0,angle,npts+1)'; % parametric angle
t(end) = [];                   % remove duplicate

tol = getfemobjectoptions('tolerancepoint');

% Build parametric cylinder
if abs(angle - 2*pi) < tol
    [X,Y,Z] = cylinder(r,npts);
    Z = Z * h;
    nodecoord_bot = [X(1,:)', Y(1,:)', Z(1,:)'];
    nodecoord_top = [X(2,:)', Y(2,:)', Z(2,:)'];
else
    x = r * cos(t);
    y = r * sin(t);
    z_bot = zeros(size(t));
    z_top = h * ones(size(t));
    nodecoord_bot = [x, y, z_bot];
    nodecoord_top = [x, y, z_top];
end

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
nodecoord_bot = nodecoord_bot * R + c;
nodecoord_top = nodecoord_top * R + c;
center_base = [0, 0, 0] * R + c;
center_top  = [0, 0, h] * R + c;

% Plot using patch
options = patchoptions(C.indim,varargin{:});
holdState = ishold;
hold on
H = struct();

% Bottom and top rings
H.rings = gobjects(2,1);
if abs(angle - 2*pi) < tol
    % Full circle: use patch for closed loop
    connec = [1:npts,1]; % connectivity for a closed loop
    H.rings(1) = patch('Faces',connec,'Vertices',nodecoord_bot,options{:}); % bottom circle
    H.rings(2) = patch('Faces',connec,'Vertices',nodecoord_top,options{:}); % top circle
else
    % Partial cylinder: create a sector (fan) for open circle arc by including center point
    nodecoord_bot_fan = [center_base; nodecoord_bot]; % add base center
    nodecoord_top_fan = [center_top; nodecoord_top]; % add top center
    connec_fan = 1:(npts+1);
    H.rings(1) = patch('Faces',connec_fan,'Vertices',nodecoord_bot_fan,options{:});
    H.rings(2) = patch('Faces',connec_fan,'Vertices',nodecoord_top_fan,options{:});
end

% Vertical spoke along the central axis (center base to center top)
H.axis = patch('Faces',[1 2],'Vertices',[center_base; center_top],options{:},'LineStyle','-.');

% Vertical spoke lines and radial lines
spoke_angles = [0, pi/2, pi, 3*pi/2];
spoke_idx = 1; % always include start
for k=2:4
    if angle >= spoke_angles(k)-tol
        [~,idx] = min(abs(t - spoke_angles(k)));
        spoke_idx = [spoke_idx, idx];
    end
end
if abs(angle - 2*pi) > tol
    spoke_idx = unique([spoke_idx, length(t)]); % also include end for open arc
end
spoke_idx = unique(spoke_idx);

H.spokes  = gobjects(length(spoke_idx),1);
H.radials = gobjects(length(spoke_idx),2);
for k = 1:length(spoke_idx)
    idx = spoke_idx(k);
    H.spokes(k) = patch('Faces',[1 2],'Vertices',[nodecoord_bot(idx,:); nodecoord_top(idx,:)],options{:});
    H.radials(k,1) = patch('Faces',[1 2],'Vertices',[center_base; nodecoord_bot(idx,:)],options{:},'LineStyle',':');
    H.radials(k,2) = patch('Faces',[1 2],'Vertices',[center_top; nodecoord_top(idx,:)],options{:},'LineStyle',':');
end
if ~holdState
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

if nargout>=1
    varargout{1} = H;
end
