function varargout = plot(C,varargin)
% function varargout = plot(C,varargin)

npts = getcharin('npts',varargin,200);
t = linspace(0,2*pi,npts);
t = t(1:end);

x = C.a*cos(t');
y = C.b*sin(t');

switch C.indim
    case 2
        nodecoord = [x,y];
        
        v = [C.vx,C.vy];
        v = v/norm(v);
        R = [v(1) v(2);
            -v(2) v(1)];
        
        c = [C.cx,C.cy];
    case 3
        nodecoord = [x,y,zeros(length(t),1)];
        
        n = [C.nx,C.ny,C.nz];
        n = n/norm(n);
        Q = [0 -n(3) n(2);
            n(3) 0 -n(1);
            -n(2) n(1) 0];
        v = [C.vx,C.vy];
        v = v/norm(v);
        R = eye(3)-v(2)*Q+(1-v(1))*Q^2;
        
        c = [C.cx,C.cy,C.cz];
end

nodecoord = nodecoord*R + repmat(c,length(t),1);
connec = [1:length(t),1];
options = patchoptions(C.indim,varargin{:});

H = patch('faces',connec,'vertices',nodecoord,options{:});

if ~(C.indim==3 && all(nodecoord(:,3)==C.cz))
    axis image
end

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
