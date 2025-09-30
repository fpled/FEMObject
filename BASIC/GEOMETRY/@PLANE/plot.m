function varargout = plot(D,varargin)
% function varargout = plot(D,varargin)

P1 = double(getcoord(D.P{1}));
P2 = double(getcoord(D.P{2}));
P3 = double(getcoord(D.P{3}));
nodecoord = [P1;P2;P3] ;
connec = [1,2,3];
edgecolor = getcharin('Color',varargin,'k') ;

% Plot using patch
hs = ishold;
hold on

H = patch('Faces',connec,'Vertices',nodecoord,'EdgeColor',edgecolor);

if ~hs
    hold off
end

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

if nargout
    varargout{1} = H;
end
