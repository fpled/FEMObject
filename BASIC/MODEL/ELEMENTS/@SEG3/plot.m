function varargout = plot(elem,node,varargin)
% function varargout = plot(elem,node,varargin)

nodecoord = double(getcoord(node));
[a,connec]=ismember(getconnec(elem),getnumber(node));

n = size(nodecoord);
if n(2)==1
    nodecoord = [nodecoord,zeros(n(1),1)];
end

options = varargin;

colelem = getcharin('facevertexcdata',options);
facecolor = getcharin('facecolor',options,'none');
edgecolor = getcharin('edgecolor',options,'k');
if ~strcmpi(facecolor,'none')
    edgecolor = facecolor;
end

if size(colelem,1)==1 && ~strcmpi(facecolor,'flat')
    colelem = colelem*ones(size(nodecoord,1),1);
elseif strcmpi(facecolor,'flat')
    nodecoord = [nodecoord(connec(:,1),:);nodecoord(connec(:,2),:);...
        nodecoord(connec(:,3),:)];
    connec = [1:size(connec,1);size(connec,1)+1:2*size(connec,1);...
        2*size(connec,1)+1:3*size(connec,1)]';
    if size(colelem,1)==1
        colelem = colelem*ones(size(nodecoord,1),1);
    else
        colelem = [colelem;colelem];
    end
end

options = setcharin('edgecolor',options,edgecolor);
if ~isempty(colelem)
    options = setcharin('facevertexcdata',options,colelem);
end

H(1) = patch('faces',connec(:,1:2),'vertices',nodecoord,options{:});
H(2) = patch('faces',connec(:,2:3),'vertices',nodecoord,options{:});

if nargout>=1
    varargout{1} = H;
end
