function varargout = plot(D,varargin)
% function varargout = plot(D,varargin)

P = vertcat(D.P{:});
nodecoord = permute(double(getcoord(P)),[3,2,1]);
connec = 1:4;
options = patchoptions(3,varargin{:});

H = patch('faces',connec,'vertices',nodecoord,options{:});

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

if nargout>=1
    varargout{1} = H;
end
