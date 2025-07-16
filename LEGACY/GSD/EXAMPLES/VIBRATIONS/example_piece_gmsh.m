%% importation du modele
S = gmsh2femobject(3,'C:\PROGRAMMES\gmsh\demos\piece.msh',3);

%% facets 
facegroup{1} = 1; % (surface interne - bas axe)
facegroup{2} = 2; % (face du dessous)
facegroup{3} = 3:6; % (surface interne - exterieur axe)
facegroup{4} = 7:10; % (surface interieure milieu axe)
facegroup{5} = 11 ; % (surface interne - haut axe)
facegroup{6} = 12:46 ; % (contour exterieur)
facegroup{7} = 47; % (face du dessus)
facegroup{8} = 48:51; % (contour exterieur du haut de l'axe)
facegroup{9} = 52:55; % (contour interieur du haut de l'axe)
facegroup{10} = 56; % (surface superieure de l'axe)
facegroup{11} = 57:60; % (contour exterieur du bas de l'axe)
facegroup{12} = 61:64; % (contour interieur du bas de l'axe)
facegroup{13} = 65; % (surface inferieure de l'axe)
facesinternes = [facegroup{[4,9,12]}];
facesefforts = [29 30 31 32];
% facesefforts = [26,31];
facefict = [facegroup{[1,3,5]}];
% facesdirichlet = [facesexternes,facesinternes];
facesdirichlet = [facesinternes];

%% construction du modele EF et CL
mat = MATERIALS(ELAS_ISOT('E',1,'NU',0.3,'RHO',1));
S = setmaterial(S,mat{1});
S = splitgroupelem(S,3000);
S = final(S);

nodesdirichlet=[];
for l=facesdirichlet % 20:23
    nodesdirichlet = union(nodesdirichlet,getnumber(getnode(getfacet(S,l))));
end

nodesefforts=[];
for l=facesefforts % 40
    nodesefforts = union(nodesefforts,getnumber(getnode(getfacet(S,l))));
end

S = addcl(S,nodesdirichlet,'U');

f = nodalload(S,nodesefforts,'FZ',1);

K1 = calc_rigi(S);
M1 = calc_mass(S);

%% output
numpeaks = [25,34];
nodeoutput = [];
for l=numpeaks
    nodeoutput = union(nodeoutput,getnumber(getnode(getpeak(S,l))));
end
ddloutput = findddl(S,'all',nodeoutput,'free');

%% modele amortissement proportionnel
alpha = @(w) 0.001;
beta = @(w) 0.001;

%% affichage modele
figure(1)
clf
plotfacets(S,'EdgeColor','none','FaceColor','b','FaceAlpha',1,'legend',false)
light('Position',[3 0 3])

figure(2)
clf
plotridges(S,'EdgeColor','k','FaceColor','k','EdgeAlpha',1,'FaceAlpha',1,'legend',false)
h1 = plotfacets(S,facesdirichlet,'EdgeColor','none','FaceColor','r','FaceAlpha',1);
h2 = plotfacets(S,facesefforts,'EdgeColor','none','FaceColor','y','FaceAlpha',1);
h3 = plotpeaks(S,numpeaks,'MarkerSize',20);
legend([h1(1),h2(2),h3(1)],'Dirichlet','Neumann','Output')

%% calcul vibrations sans amort
scanm = 1:30
[V,s] = calc_mode(K1,M1,scanm);
s = diag(s);
wres = sqrt(s)

mm = 7;
wband = [wres(1)/10,wres(mm)+(wres(mm+1)-wres(mm))/2];

%% affichage modes
V = unfreevector(S,V);
Vx = V(findddl(S,'UX'),:);
Vy = V(findddl(S,'UY'),:);

figure(10)
clf
scanmode = 1:6;
nbm = length(scanmode); 
for k=scanmode
    subplot(ceil(sqrt(nbm)),ceil(nbm/ceil(sqrt(nbm))),find(k==scanmode))
    title(['w_' num2str(k) ' = '  num2str(wres(k))]);
    scal = 1/max(max(abs(V(:,k))))*max(max(abs(getcoord(getnode(S)))))*0.1;
    plot(FENODEFIELD(sqrt(Vx(:,k).^2+Vy(:,k).^2)),S+scal*V(:,k));
    pause(0.5)
end
