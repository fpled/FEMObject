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
        
        % Rotation matrix
        v = [E.vx, E.vy];
        R = calcrotation(E,v);
        
        % Center
        c = [E.cx, E.cy];
        
    case 3
        z = zeros(npts,1);
        nodecoord = [x, y, z];
        
        % Rotation matrix
        v = [E.vx, E.vy];
        n = [E.nx, E.ny, E.nz];
        R = calcrotation(E,v,n);
        
        % Center
        c = [E.cx, E.cy, E.cz];
        
    otherwise
        error('Wrong space dimension');
end

% Rotate into global frame and translate to center
nodecoord = nodecoord * R + c;

% Connectivity for a closed loop
connec = [1:npts,1];

% Plot using patch
options = patchoptions(E.indim,varargin{:});
hs = ishold;
hold on

H = patch('Faces',connec,'Vertices',nodecoord,options{:});

if ~hs
    hold off
end

tol = getfemobjectoptions('tolerancepoint');
if ~(E.indim==3 && all(abs(nodecoord(:,3) - c(3)) < tol))
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

if nargout
    varargout{1} = H;
end
