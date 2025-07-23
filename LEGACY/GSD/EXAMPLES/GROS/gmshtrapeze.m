function varargout = gmshtrapeze(dim,lc1,lc2,varargin)
% function varargout = gmshtrapeze(dim,lc1,lc2,varargin)

G = GMSHFILE();
P{1} = [0,0,0];
P{2} = [2,.3,0];
P{3} = [2,.7,0];
P{4} = [0,1,0];
G = createpoints(G,P(1:4),lc1,1:4);
G = createcontour(G,1:4,1:4,100);

r=.25;
cx=.5;
cy=.5;
P{5} = [cx,cy,0];
P{6} = [cx-r,cy,0];
P{7} = [cx,cy-r,0];
P{8} = [cx+r,cy,0];
P{9} = [cx,cy+r,0];

G = createpoints(G,P(5:9),lc2,5:9);
G = createcirclecontour(G,5,6:9,5:8,101);

r=.15;
cx=1.3;
cy=.5;
P{10} = [cx,cy,0];
P{11} = [cx-r,cy,0];
P{12} = [cx,cy-r,0];
P{13} = [cx+r,cy,0];
P{14} = [cx,cy+r,0];

G = createpoints(G,P(10:14),lc2,10:14);
G = createcirclecontour(G,10,11:14,9:12,102);

G = createcurveloop(G,[1:4,-(5:8),-(9:12)],200);

G = createplanesurface(G,200,1);

n = max(nargout,1);
varargout = cell(1,n);
[varargout{:}] = gmsh2femobject(dim,G,dim:-1:dim-n+1,varargin{:});


