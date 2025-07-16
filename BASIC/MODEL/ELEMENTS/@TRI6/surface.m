function varargout = surface(elem,node,varargin)
% function varargout = surface(elem,node,varargin)

nodecoord = double(getcoord(node));
colnode = getcharin('FaceVertexCData',varargin);
nodecoord = [nodecoord,colnode];
connec = calc_conneclocal(elem,node);
connec = connec(:,[1,4,2,5,3,6]);
% connec = connec(:,[1,2,3]);

H = patch('Faces',connec,'Vertices',nodecoord,'FaceLighting','phong',varargin{:});

if nargout>=1
    varargout{1} = H;
end
