G = GMSHFILE();
G = setfile(G,'gmsh_multi_inclusions_8');

cl = .02;
D = DOMAIN(2,[0 0],[1 1]);
G = G + gmshfile(D,cl,1:4,1:4,5);
pmax = 4;
lmax = 5;
smax = 1;
lineloop = 1:4;

% points = [1,1;1,2;1,3;2,1;2,3;3,1;3,2;3,3]/4;
points = [2,2;2,5;2,8;5,2;5,8;8,2;8,5;8,8]/10;
r = 1/8;
% r = 0.13;
C = cell(1,size(points,1));
for i=1:size(points,1);
    C{i} = CIRCLE(points(i,1),points(i,2),r);
    G = G + gmshfile(C{i},cl,pmax+1,pmax+(2:5),lmax+(2:5),lmax+1,smax+1);
    lineloop = [lineloop,-lmax-(2:5)];
    pmax = pmax+5;
    lmax = lmax+5;
    smax = smax+1;
end
rc = .1;
D0 = DOMAIN(2,1/2-[rc,rc],1/2+[rc,rc]);
G = G + gmshfile(D0,cl,pmax+(1:4),lmax+(2:5),lmax+1,smax+1);
lineloop = [lineloop,-lmax-(2:5)];
pmax = pmax+5;
lmax = lmax+5;
smax = smax+1;

G = createcurveloop(G,lineloop,lmax+1);
G = createplanesurface(G,lmax+1,1);
S = gmsh2femobject(2,G);
% figure(1)
% clf
% plotparamelem(S,'group')
