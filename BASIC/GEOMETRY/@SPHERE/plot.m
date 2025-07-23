function varargout = plot(S,varargin)
% function varargout = plot(S,varargin)

npts = getcharin('npts',varargin,100); % angular resolution
t = linspace(0,2*pi,npts+1)'; % parametric angle
t(end) = []; % remove duplicate

% Radius
r = S.r;

%% Old version
% Rotate around axis n = [nx, ny, nz] by angle of rotation phi =
% atan2(vy, vx) using tangent vector v = [vx, vy] via Rodrigues'
% rotation formula
%% New version
% Twist the XY plane about z = [0, 0, 1] by phi = atan2(vy, vx)
% using tangent vector v = [vx, vy], then tilt from z axis to normal
% vector n = [nx, ny, nz] so that the circle's normal is n regardless
% of v = [vx, vy]
v = [S.vx, S.vy];
n = [S.nx, S.ny, S.nz];
R = calcrotation(S,v,n);

% Translate to center c = [cx, cy, cz]
c = [S.cx, S.cy, S.cz];

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
options = patchoptions(S.indim, varargin{:});
holdState = ishold;
hold on

% Rings
H.rings = gobjects(size(N,1),1);
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
    ui = ui/norm(ui);
    wi = cross(ni, ui);
    
    % Parametric circle in local plane
    nodecoord = r * (cos(t)*ui + sin(t)*wi);
    
    % Rotate and translate
    nodecoord = nodecoord * R + c;
    
    % draw closed ring
    H.rings(i) = patch('Faces',connec,'Vertices',nodecoord,options{:});
end

% Main diameters (radial lines)
diam_pts = [ r  0  0;
            -r  0  0;
             0  r  0;
             0 -r  0;
             0  0  r;
             0  0 -r];

% Rotate and translate endpoints
H.diameters = gobjects(3,1);
for i = 1:3
    p1 = diam_pts(2*i-1,:) * R + c; % +axis
    p2 = diam_pts(2*i,:)   * R + c; % -axis
    H.diameters(i) = patch('Faces',[1 2],'Vertices',[p1; p2],options{:},'LineStyle',':');
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
