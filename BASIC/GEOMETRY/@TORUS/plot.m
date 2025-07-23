function varargout = plot(C,varargin)
% function varargout = plot(C,varargin)

npts = getcharin('npts',varargin,200); % angular resolution
% t = linspace(0,2*pi,npts+1)'; % parametric angle
% t(end) = []; % remove duplicate

% Radius and height
r = C.r;
h = C.h;

% Build parametric cylinder
% x = r * cos(t);
% y = r * sin(t);
% z_bot = zeros(npts,1);
% z_top = h * ones(npts,1);
% nodecoord_bot = [x, y, z_bot];
% nodecoord_top = [x, y, z_top];

% Build parametric cylinder
[X,Y,Z] = cylinder(r,npts);
Z = Z * h;
nodecoord_bot = [X(1,:)', Y(1,:)', Z(1,:)'];
nodecoord_top = [X(2,:)', Y(2,:)', Z(2,:)'];

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

% Connectivity for a closed loop
connec = [1:npts,1];

% Plot using patch
options = patchoptions(C.indim,varargin{:});
holdState = ishold;
hold on
% bottom ring
Hb = patch('Faces',connec,'Vertices',nodecoord_bot,options{:});
% top ring
Ht = patch('Faces',connec,'Vertices',nodecoord_top,options{:});
% four vertical spokes at t = 0, pi/2, pi, 3*pi/2
idx = round(linspace(1,npts,5));  idx(end)=[];  % [0, pi/2, pi, 3*pi/2]
Hs = gobjects(4,1);
for k=1:4
    nodecoord = [nodecoord_bot(idx(k),:);
                 nodecoord_top(idx(k),:)];
    Hs(k) = patch('Faces',[1 2],'Vertices',nodecoord,options{:});
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
