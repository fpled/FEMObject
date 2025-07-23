function varargout = plot(E,varargin)
% function varargout = plot(E,varargin)

npts = getcharin('npts',varargin,200); % angular resolution
t = linspace(0,2*pi,npts+1)'; % parametric angle
t(end) = []; % remove duplicate

% Semi-axes
a = E.a;
b = E.b;

% Build parametric ellipse
x = a * cos(t);
y = b * sin(t);

switch E.indim
    case 2
        nodecoord = [x, y];
        
        % Rotate in-plane around z = [0, 0, 1] by angle of rotation theta =
        % atan2(vy, vx) using tangent vector v = [vx, vy]
        v = [E.vx, E.vy];
        R = calcrotation(E,v);
        
        % Translate to center c = [cx, cy]
        c = [E.cx, E.cy];
        
    case 3
        nodecoord = [x, y, zeros(npts,1)];
        
        %% Old version
        % Rotate around axis n = [nx, ny, nz] by angle of rotation phi =
        % atan2(vy, vx) using tangent vector v = [vx, vy] via Rodrigues'
        % rotation formula
        %% New version
        % Twist the XY plane about z = [0, 0, 1] by phi = atan2(vy, vx)
        % using tangent vector v = [vx, vy], then tilt from z axis to normal
        % vector n = [nx, ny, nz] so that the circle's normal is n regardless
        % of v = [vx, vy]
        v = [E.vx, E.vy];
        n = [E.nx, E.ny, E.nz];
        R = calcrotation(E,v,n);
        
        % Translate to center c = [cx, cy, cz]
        c = [E.cx, E.cy, E.cz];
        
    otherwise
        error('Wrong space dimension');
end

% Rotate and translate
nodecoord = nodecoord * R + c;

% Connectivity for a closed loop
connec = [1:npts,1];

% Plot using patch
options = patchoptions(E.indim,varargin{:});
H = patch('Faces',connec,'Vertices',nodecoord,options{:});

tol = getfemobjectoptions('tolerancepoint');
if ~(E.indim==3 && all(abs(nodecoord(:,3) - E.cz) < tol))
    axis image
end

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

if nargout>=1
    varargout{1} = H;
end
