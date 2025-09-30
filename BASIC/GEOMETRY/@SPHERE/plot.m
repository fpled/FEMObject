function varargout = plot(S,varargin)
% function varargout = plot(S,varargin)

npts = getcharin('npts',varargin,200); % angular resolution
t = linspace(0,2*pi,npts+1)'; % parametric angle
t(end) = []; % remove duplicate

% Radius
r = S.r;

% Rotation matrix
v = [S.vx, S.vy];
n = [S.nx, S.ny, S.nz];
R = calcrotation(S,v,n);

% Center
c = [S.cx, S.cy, S.cz];

% Connectivity for a closed loop
connec = [1:npts,1];

% Plot using patch
options = patchoptions(S.indim, varargin{:});
hs = ishold;
hold on

% Three principal-section circles
% First circle
x1 = r * cos(t);
y1 = r * sin(t);
z1 = zeros(npts,1);
nodecoord1 = [x1, y1, z1];

% Second circle
x2 = r * cos(t);
y2 = zeros(npts,1);
z2 = r * sin(t);
nodecoord2 = [x2, y2, z2];

% Third circle
x3 = zeros(npts,1);
y3 = r * cos(t);
z3 = r * sin(t);
nodecoord3 = [x3, y3, z3];

% Rotate into global frame and translate to center
nodecoord1 = nodecoord1 * R + c;
nodecoord2 = nodecoord2 * R + c;
nodecoord3 = nodecoord3 * R + c;

% Draw circles
H.circles = gobjects(3,1);
H.circles(1) = patch('Faces',connec,'Vertices',nodecoord1,options{:});
H.circles(2) = patch('Faces',connec,'Vertices',nodecoord2,options{:});
H.circles(3) = patch('Faces',connec,'Vertices',nodecoord3,options{:});

% Diameters
diamcoord = [ r  0  0;
             -r  0  0;
              0  r  0;
              0 -r  0;
              0  0  r;
              0  0 -r];

% Rotate into global frame and translate to center
diamcoord = diamcoord * R + c;

% Draw diameters
H.diameters = gobjects(3,1);
for i = 1:3
    p1 = diamcoord(2*i-1,:);
    p2 = diamcoord(2*i,:);
    H.diameters(i) = patch('Faces',[1 2],'Vertices',[p1; p2],options{:},'LineStyle',':');
end

H = [H.circles(:); H.diameters(:)];

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
