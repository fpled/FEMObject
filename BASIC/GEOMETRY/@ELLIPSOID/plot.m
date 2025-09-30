function varargout = plot(E,varargin)
% function varargout = plot(E,varargin)

npts = getcharin('npts',varargin,200); % angular resolution
t = linspace(0,2*pi,npts+1)'; % parametric angle
t(end) = []; % remove duplicate

% Semi-axes
a = E.a;
b = E.b;
c = E.c;

% Rotation matrix
v = [E.vx, E.vy];
n = [E.nx, E.ny, E.nz];
R = calcrotation(E,v,n);

% Center
ct = [E.cx, E.cy, E.cz];

% Connectivity for a closed loop
connec = [1:npts,1];

% Plot using patch
options = patchoptions(E.indim, varargin{:});
hs = ishold;
hold on

H = struct();

% Three principal-section ellipses
% First ellipse
x1 = a * cos(t);
y1 = b * sin(t);
z1 = zeros(npts,1);
nodecoord1 = [x1, y1, z1];

% Second ellipse
x2 = a * cos(t);
y2 = zeros(npts,1);
z2 = c * sin(t);
nodecoord2 = [x2, y2, z2];

% Third ellipse
x3 = zeros(npts,1);
y3 = b * cos(t);
z3 = c * sin(t);
nodecoord3 = [x3, y3, z3];

% Rotate into global frame and translate to center
nodecoord1 = nodecoord1 * R + ct;
nodecoord2 = nodecoord2 * R + ct;
nodecoord3 = nodecoord3 * R + ct;

% Draw ellipses
H.ellipses = gobjects(3,1);
H.ellipses(1) = patch('Faces',connec,'Vertices',nodecoord1,options{:});
H.ellipses(2) = patch('Faces',connec,'Vertices',nodecoord2,options{:});
H.ellipses(3) = patch('Faces',connec,'Vertices',nodecoord3,options{:});

% Diameters
diamcoord = [ a  0  0;
             -a  0  0;
              0  b  0;
              0 -b  0;
              0  0  c;
              0  0 -c];

% Rotate into global frame and translate to center
diamcoord = diamcoord * R + ct;

H.diameters = gobjects(3,1);
for i = 1:3
    p1 = diamcoord(2*i-1,:);
    p2 = diamcoord(2*i,:);
    % Draw diameters
    H.diameters(i) = patch('Faces',[1 2],'Vertices',[p1; p2],options{:},'LineStyle',':');
end

H = [H.ellipses(:); H.diameters(:)];

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
elseif E.indim==3
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
