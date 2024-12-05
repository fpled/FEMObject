function varargout = plot(D,varargin)
% function varargout = plot(D,varargin)

P = double(getcoord(D.P));
nodecoord = reshape(P,getindim(D.P),numel(D.P))';
nodecoord = [nodecoord;nodecoord(end,:)];
connec = 1:size(nodecoord,1);
options = patchoptions(2,'facevertexcdata',1,varargin{:});

H = patch('faces',connec,'vertices',nodecoord,options{:});

if nargout>=1
    varargout{1} = H;
end
  