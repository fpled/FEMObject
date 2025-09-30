function varargout = plot(D,varargin)
% function varargout = plot(D,varargin)

P = double(getcoord(D.P));
nodecoord = reshape(P,getindim(D.P),numel(D.P))';
nodecoord = [nodecoord;nodecoord(end,:)];
connec = 1:size(nodecoord,1);

% Plot using patch
options = patchoptions(2,'FaceVertexCData',1,varargin{:});
hs = ishold;
hold on

H = patch('Faces',connec,'Vertices',nodecoord,options{:});

if ~hs
    hold off
end

% Optional view or camera controls
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

if nargout
    varargout{1} = H;
end
  