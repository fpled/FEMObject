function G = gmshBeam(P,cl,filename,indim,varargin)
% function G = gmshBeam(P,cl,filename,indim,varargin)
% P : cell array
% cl : characteristic length
% filename : file name (optional)
% indim : space dimension (optional, max(cellfun(@(x) size(x,2),P)) by default)

if nargin<4 || isempty(indim)
    indim = max(cellfun(@(x) size(x,2),P));
end

dim = 1;

G = GMSHFILE();
if nargin>=3 && ischar(filename)
    G = setfile(G,filename);
end

numberpoints = 1:numel(P);
numberlines = numberpoints(1:end-1);
seg = [1:length(numberpoints)-1;2:length(numberpoints)]';
seg = numberpoints(seg);
G = createpoints(G,P,cl,numberpoints);
G = createlines(G,seg,numberlines);

G = gmsh2femobject(indim,G,dim,varargin{:});
