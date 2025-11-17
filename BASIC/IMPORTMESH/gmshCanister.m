function varargout = gmshCanister(cl1,cl2,cl0,cltip,filename,indim,varargin)
% function varargout = gmshCanister(cl1,cl2,cl0,cltip,filename,indim)
% cl1, cl2, cl0, cltip : characteristic lengths
% filename : file name (optional)
% indim : space dimension (optional, 2 by default)

if nargin<6 || isempty(indim)
    indim = 2;
end

dim = 2;

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
if nargin>=5 && ischar(filename)
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

G = createlines(G,[[13 15];[15 16];[16 14]],16:18);
G = createcurveloop(G,[-14,16:18],3);

G = createplanesurface(G,1,1);
G = createplanesurface(G,2,2);
G = createplanesurface(G,3,3);

if ischarin('recombine',varargin)
    G = recombinesurface(G);
end
varargin = delonlycharin('recombine',varargin);

n = max(nargout,1);
varargout = cell(1,n);
[varargout{:}] = gmsh2femobject(indim,G,dim:-1:dim-n+1,varargin{:});
