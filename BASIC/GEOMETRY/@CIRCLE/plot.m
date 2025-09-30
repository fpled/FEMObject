function varargout = plot(C,varargin)
% function varargout = plot(C,varargin)

npts = getcharin('npts',varargin,200); % angular resolution
t = linspace(0,2*pi,npts+1)'; % parametric angle
t(end) = []; % remove duplicate

% Radius
r = C.r;

% Build parametric circle
x = r * cos(t);
y = r * sin(t);

switch C.indim
    case 2
        nodecoord = [x, y];
        
        % Rotation matrix
        v = [C.vx, C.vy];
        R = calcrotation(C,v);
        
        % Center
        c = [C.cx, C.cy];
        
    case 3
        z = zeros(npts,1);
        nodecoord = [x, y, z];
        
        % Rotation matrix
        v = [C.vx, C.vy];
        n = [C.nx, C.ny, C.nz];
        R = calcrotation(C,v,n);
        
        % Center
        c = [C.cx, C.cy, C.cz];
        
    otherwise
        error('Wrong space dimension');
end

% Rotate into global frame and translate to center
nodecoord = nodecoord * R + c;

% Connectivity for a closed loop
connec = [1:npts,1];

% Plot using patch
options = patchoptions(C.indim,varargin{:});
hs = ishold;
hold on

H = patch('Faces',connec,'Vertices',nodecoord,options{:});

if ~hs
    hold off
end

tol = getfemobjectoptions('tolerancepoint');
if ~(C.indim==3 && all(abs(nodecoord(:,3) - c(3)) < tol))
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

if nargout
    varargout{1} = H;
end
