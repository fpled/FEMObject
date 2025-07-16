function varargout = plot(E,varargin)
% function varargout = plot(E,varargin)

npts = getcharin('npts',varargin,200); % angular resolution
t = linspace(0,2*pi,npts+1)';          % angle
t(end) = [];                           % avoid duplicate point

% Semi-axes
a = E.a;
b = E.b;
c = E.c;

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
ct = [E.cx, E.cy, E.cz];

% Define three/six plane normals
N = [1 0 0;
     0 1 0;
     0 0 1];
% N = [1  0  0;
%      0  1  0;
%      0  0  1;
%      1  1  0;
%      0  1  1;
%      1  0  1 ];
N = N ./ vecnorm(N,2,2);  % normalize each plane normal

% Connectivity for a closed loop
connec = [1:npts,1];

% Plot using patch
options = patchoptions(E.indim, varargin{:});
holdState = ishold;
hold on

H = gobjects(size(N,1),1);
for i = 1:size(N,1)
    ni = N(i,:);

    % Pick arbitrary vector not parallel to ni
    if abs(ni(1))<0.9
        vec = [1 0 0];
    else
        vec = [0 1 0];
    end
    
    % Build local orthonormal basis (ui, wi) in plane orthogonal to ni
    ui = cross(vec, ni);
    ui = ui / norm(ui);
    wi = cross(ni, ui);
    
    % Parametric directions in local plane
    D = cos(t)*ui + sin(t)*wi;
    
    % Scale each D(k,:) so (x/a)^2+(y/b)^2+(z/c)^2 == 1
    S = 1./sqrt( (D(:,1)/a).^2 + (D(:,2)/b).^2 + (D(:,3)/c).^2 );
    nodecoord = D .* S;

    % Rotate and translate
    nodecoord = nodecoord * R + ct;
    
    % draw closed ring
    H(i) = patch('Faces',connec,'Vertices',nodecoord,options{:});
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
