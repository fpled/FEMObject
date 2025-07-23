function varargout = plot(L,varargin)
% function varargout = plot(L,varargin)

P1 = getcoord(L.P{1});
P2 = getcoord(L.P{2});
nodecoord = double([P1;P2]) ;
n = size(nodecoord);
nodecoord = [nodecoord,zeros(n(1),n(2)==1)];
connec = 1:2;
edgecolor = getcharin('Color',varargin,'k') ;

H = patch('Faces',connec,'Vertices',nodecoord,'EdgeColor',edgecolor);

if nargout>=1
    varargout{1} = H;
end
  