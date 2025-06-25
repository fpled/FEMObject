function varargout = gmshmulti(I,cl,filename,indim,varargin)
% function varargout = gmshmulti(I,cl,filename,indim)
% I : DOMAIN or CIRCLE or ELLIPSE or QUADRANGLE
% cl : characteristic length
% filename : file name (optional)
% indim : space dimension (optional, max(cellfun(@(x) getindim(x),I)) by default)

if nargin<4 || isempty(indim)
    indim = max(cellfun(@(x) getindim(x),I));
end

if ~iscell(I)
    I = {I};
end
if isscalar(cl)
    cl = repmat(cl,1,length(I));
end

G = GMSHFILE();
if nargin>=3 && ischar(filename)
    G = setfile(G,filename);
end

numpoints = 1:5;
numlines = 1:5;
% numlineloop = 1:4;
for j=1:length(I)
%     numlineloop = [numlineloop,-numlines(1:end-1)];
    if isa(I{j},'DOMAIN') || isa(I{j},'QUADRANGLE')
        GI = gmshfile(I{j},cl(j),numpoints(1:end-1),numlines(1:end-1),numlines(end),j+1,varargin{:});
    elseif isa(I{j},'CIRCLE') || isa(I{j},'ELLIPSE')
        GI = gmshfile(I{j},cl(j),numpoints(1),numpoints(2:end),numlines(1:end-1),numlines(end),j+1,varargin{:});
    end
    G = G+GI;
    numpoints = numpoints+5;
    numlines = numlines+5;
end
% G = createcurveloop(G,numlineloop,numlines(end));
% G = createplanesurface(G,numlines(end),1);
% if ischarin('recombine',varargin)
%     G = recombinesurface(G,1);
% end
varargin = delonlycharin('recombine',varargin);

n=max(nargout,1);
varargout = cell(1,n);
dim = max(cellfun(@(x) getdim(x),I));
[varargout{:}] = gmsh2femobject(indim,G,dim:-1:dim-n+1,varargin{:});
