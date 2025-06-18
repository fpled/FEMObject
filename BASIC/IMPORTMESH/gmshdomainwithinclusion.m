function varargout = gmshdomainwithinclusion(D,I,clD,clI,filename,indim,varargin)
% function varargout = gmshdomainwithinclusion(D,I,clD,clI,filename,indim)
% D : DOMAIN or QUADRANGLE
% I : DOMAIN or CIRCLE or ELLIPSE or QUADRANGLE or LIGNE or POINT
% clD, clI : characteristic lengths
% filename : file name (optional)
% indim : space dimension (optional, getindim(D) by default)

if nargin<6 || isempty(indim)
    indim = getindim(D);
end
if nargin<4 || isempty(clI)
    clI = clD;
end

if ~iscell(I)
    I = {I};
end
if length(clI)==1
    clI = repmat(clI,1,length(I));
end

G = gmshfile(D,clD,1:4,1:4,5);
if nargin>=5 && ischar(filename)
    G = setfile(G,filename);
end

numpoints = 4+(1:5);
numlines = 5+(1:5);
numlineloop = 1:4;
numberembeddedpoints = [];
numberembeddedlines = [];
for j=1:length(I)
    if isa(I{j},'POINT')
        GI = gmshfile(I{j},clI(j),numpoints(1));
        numberembeddedpoints = [numberembeddedpoints,numpoints(1)];
    elseif isa(I{j},'LIGNE')
        GI = gmshfile(I{j},clI(j),numpoints(1:2),numlines(1));
        numberembeddedlines = [numberembeddedlines,numlines(1)];
    else
        if isa(I{j},'DOMAIN') || isa(I{j},'QUADRANGLE')
            GI = gmshfile(I{j},clI(j),numpoints(1:end-1),numlines(1:end-1),numlines(end),j+1,varargin{:});
        elseif isa(I{j},'CIRCLE') || isa(I{j},'ELLIPSE')
            GI = gmshfile(I{j},clI(j),numpoints(1),numpoints(2:end),numlines(1:end-1),numlines(end),j+1,varargin{:});
        end
        GI = createphysicalsurface(GI,j+1,j+1);
        numlineloop = [numlineloop,-numlines(1:end-1)];
    end
    G = G+GI;
    numpoints = numpoints+5;
    numlines = numlines+5;
end
G = createcurveloop(G,numlineloop,numlines(end));
G = createplanesurface(G,numlines(end),1);
if ischarin('recombine',varargin)
    G = recombinesurface(G,1);
end
if ~isempty(numberembeddedpoints)
    G = embedpointsinsurface(G,numberembeddedpoints,1);
end
if ~isempty(numberembeddedlines)
    G = embedcurvesinsurface(G,numberembeddedlines,1);
end
G = createphysicalsurface(G,1,1);
varargin = delonlycharin('recombine',varargin);

n=max(nargout,1);
varargout = cell(1,n);
[varargout{:}] = gmsh2femobject(indim,G,getdim(D):-1:getdim(D)-n+1,varargin{:});
