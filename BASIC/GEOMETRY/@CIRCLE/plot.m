function varargout = plot(C,varargin)
% function varargout = plot(C,varargin)

npts = getcharin('npts',varargin,200); % angular resolution
t = linspace(0,2*pi,npts+1)';          % angle
t(end) = [];                           % avoid duplicate point

% Radius
r = C.r;

% Build parametric circle
x = r * cos(t);
y = r * sin(t);

switch C.indim
    case 2
        nodecoord = [x, y];
        
        % Rotate in-plane around z = [0, 0, 1] by angle of rotation theta =
        % atan2(vy, vx) using tangent vector v = [vx, vy]
        v = [C.vx, C.vy];
        R = calcrotation(C,v);
        
        % Translate to center c = [cx, cy]
        c = [C.cx, C.cy];
        
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
        v = [C.vx, C.vy];
        n = [C.nx, C.ny, C.nz];
        R = calcrotation(C,v,n);
        
        % Translate to center c = [cx, cy, cz]
        c = [C.cx, C.cy, C.cz];
        
    otherwise
        error('Wrong space dimension');
end

% Rotate and translate
nodecoord = nodecoord * R + c;

% Connectivity for a closed loop
connec = [1:npts,1];

% Plot using patch
options = patchoptions(C.indim,varargin{:});
H = patch('Faces',connec,'Vertices',nodecoord,options{:});

tol = getfemobjectoptions('tolerancepoint');
if ~(C.indim==3 && all(abs(nodecoord(:,3) - C.cz) < tol))
    axis image
end

% Optional view or camera controls
numview = getcharin('view',varargin);
up_vector = getcharin('camup',varargin);
camera_position = getcharin('campos',varargin);
if ~isempty(numview)
    view(numview)
elseif C.indim==3
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
