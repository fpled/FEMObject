function varargout = plot(D,varargin)
% function varargout = plot(D,varargin)

P = vertcat(D.P{:});
nodecoord = permute(double(getcoord(P)),[3,2,1]);
connec = 1:4;

% Plot using patch
options = patchoptions(3,varargin{:});
hs = ishold;
hold on

H = patch('Faces',connec,'Vertices',nodecoord,options{:});

if ~hs
    hold off
end

axis image

numview = getcharin('view',varargin);
up_vector = getcharin('camup',varargin);
camera_position = getcharin('campos',varargin);
if ~isempty(numview)
    view(numview)
elseif D.indim==3
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
