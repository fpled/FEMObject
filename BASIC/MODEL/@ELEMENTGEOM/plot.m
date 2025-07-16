function varargout = plot(elem,node,varargin)
% function varargout = plot(elem,node,varargin)

if nargin==2
    options = patchoptions(getdim(elem));
elseif getdim(elem)==0
    options = patchoptions(getdim(elem),varargin{:});
else
    options = varargin;
end

nodecoord = double(getcoord(node));
connec = calc_conneclocal(elem,node);

faces = patchfaces(elem,connec);
facevertexcdata = getcharin('FaceVertexCData',options);

if ~isempty(facevertexcdata) && size(facevertexcdata,1)==getnbelem(elem)
    facevertexcdata = repmat(facevertexcdata',size(faces,1)/getnbelem(elem),1);
    facevertexcdata = facevertexcdata(:);
    options = setcharin('FaceVertexCData',options,facevertexcdata);
end

H = patch('Faces',faces,'Vertices',nodecoord,options{:});

if nargout>=1
    varargout{1} = H;
end
