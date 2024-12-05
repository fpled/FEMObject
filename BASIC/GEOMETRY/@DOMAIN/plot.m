function varargout = plot(D,varargin)
% function varargout = plot(D,varargin)

P1 = double(getcoord(D.P1));
P2 = double(getcoord(D.P2));

switch D.indim
    case 1
        nodecoord = [P1;P2] ;
        n = size(nodecoord);
        nodecoord = [nodecoord,zeros(n(1),n(2)==1)];
        connec = 1:2;
        
    case 2
        nodecoord = [P1(1),P1(2);P2(1),P1(2);P2(1),P2(2);P1(1),P2(2)];
        connec = 1:4;
        
    case 3
        nodecoord = [P1(1),P1(2),P1(3);P2(1),P1(2),P1(3);...
            P2(1),P2(2),P1(3);P1(1),P2(2),P1(3);...
            P1(1),P1(2),P2(3);P2(1),P1(2),P2(3);...
            P2(1),P2(2),P2(3);P1(1),P2(2),P2(3)];
        connec = [1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 1 2 3 4; 5 6 7 8];

    otherwise
        error('Wrong space dimension')
        
end

options = patchoptions(D.indim,'facevertexcdata',1,varargin{:});

H = patch('faces',connec,'vertices',nodecoord,options{:});

axis image

numview = getcharin('view',varargin);
up_vector = getcharin('camup',varargin);
if ~isempty(numview)
    view(numview)
elseif D.indim==3
    view(3)
end
if ~isempty(up_vector)
    camup(up_vector)
end

if nargout>=1
    varargout{1} = H;
end
