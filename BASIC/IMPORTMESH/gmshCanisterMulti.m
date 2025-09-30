function varargout = gmshCanisterMulti(I,cl1,cl2,cl0,cltip,clI,filename,indim,varargin)
% function varargout = gmshCanisterMulti(I,cl1,cl2,cl0,cltip,clI,filename,indim)
% I : DOMAIN or CIRCLE or ELLIPSE or QUADRANGLE
% cl1, cl2, cl0, cltip, clI : characteristic lengths
% filename : file name (optional)
% indim : space dimension (optional, 2 by default)

if nargin<8 || isempty(indim)
    indim = 2;
end
if nargin<6 || isempty(clI)
    clI = cl1;
end

if ~iscell(I)
    I = {I};
end
if isscalar(clI)
    clI = repmat(clI,1,length(I));
end

lin = 0.3;
lmid = 0.5;
lout = 0.3;
t = 0.02;
hin = 0.15;
hout = 0.3;
h1 = 1;
h0 = 0.1;
h2 = 0.5;
P{1} =  [0            , 0       , 0];
P{2} =  [lout+lmid-t  , 0       , 0];
P{3} =  [lout+lmid-t  , h1-hin  , 0];
P{4} =  [lout+lmid    , h1-hin  , 0];
P{5} =  [lout+lmid    , 0       , 0];
P{6} =  [lout+lmid+lin, 0       , 0];
P{7} =  [lout+lmid+lin, h1      , 0];
P{8} =  [lout+t       , h1      , 0];
P{9} =  [lout+t       , hout    , 0];
P{10} = [lout         , hout    , 0];
P{11} = [lout         , h1      , 0];
P{12} = [0            , h1      , 0];
P{13} = [lout         , h1+h0   , 0];
P{14} = [0            , h1+h0   , 0];
P{15} = [lout         , h1+h0+h2, 0];
P{16} = [0            , h1+h0+h2, 0];

G = GMSHFILE();
if nargin>=7 && ischar(filename)
    G = setfile(G,filename);
end

G = createpoints(G,P([1 2 7 8]),cl2,[1 2 7 8]);
G = createpoints(G,P([3 4 9 10]),cltip,[3 4 9 10]);
G = createpoints(G,P(5:6),cl1,5:6);
G = createpoints(G,P(11:14),cl0,11:14);
G = createpoints(G,P(15:16),cl2,15:16);
G = createcontour(G,1:12,1:12,1);

G = createlines(G,[[11 13];[13 14];[14 12]],13:15);
G = createcurveloop(G,[-11,13:15],2);
G = createplanesurface(G,2,2);
if ischarin('recombine',varargin)
    G = recombinesurface(G,2);
end

G = createlines(G,[[13 15];[15 16];[16 14]],16:18);
G = createcurveloop(G,[-14,16:18],3);
G = createplanesurface(G,3,3);
if ischarin('recombine',varargin)
    G = recombinesurface(G,3);
end

numpoints = 16;
numlines = 18;
numlineloop = 3;
numsurface = 3;
numcurves = 1:12;
for j=1:length(I)
    if isa(I{j},'DOMAIN') || isa(I{j},'QUADRANGLE')
        numpoints = numpoints(end)+(1:4);
        numlines = numlines(end)+(1:4);
        numlineloop = numlineloop(end)+1;
        numsurface = numsurface(end)+1;
        GI = gmshfile(I{j},clI(j),numpoints,numlines,numlineloop,numsurface,varargin{:});
    elseif isa(I{j},'CIRCLE') || isa(I{j},'ELLIPSE')
        numpoints = numpoints(end)+(1:5);
        numlines = numlines(end)+(1:4);
        numlineloop = numlineloop(end)+1;
        numsurface = numsurface(end)+1;
        GI = gmshfile(I{j},clI(j),numpoints(1),numpoints(2:end),numlines,numlineloop,numsurface,varargin{:});
    end
    numcurves = [numcurves,-numlines];
    G = G+GI;
end

numlineloop = numlineloop(end)+1;
G = createcurveloop(G,numcurves,numlineloop);
G = createplanesurface(G,numlineloop,1);
if ischarin('recombine',varargin)
    G = recombinesurface(G,1);
end

varargin = delonlycharin('recombine',varargin);

n = max(nargout,1);
varargout = cell(1,n);
[varargout{:}] = gmsh2femobject(indim,G,2:-1:2-n+1,varargin{:});
