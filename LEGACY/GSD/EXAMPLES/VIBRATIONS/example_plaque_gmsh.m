%% importation du modele
Dom = QUADRANGLE([0,0,0],[1,0,0],[1,0.2,0],[0,0.4,0]);
S = gmsh(Dom,0.03);

%% construction du modele EF et CL
facesdirichlet = 4;
facesefforts = 2;

mat = MATERIALS(ELAS_SHELL('E',1,'NU',0.3,'DIM3',1,'RHO',1));
S = setmaterial(S,mat{1});
S = convertelem(S,'DKT');
S = final(S);
S = addcl(S,getedge(Dom,facesdirichlet),'U');
S = addcl(S,getedge(Dom,facesdirichlet),'R');

f = surfload(S,getedge(Dom,facesefforts),'FZ',1);

K1 = calc_rigi(S);
M1 = calc_mass(S);

%% output
numridges = [3];
nodeoutput = [];
for l=numridges
    nodeoutput = union(nodeoutput,getnumber(getnode(getridge(S,l))));
end
ddloutput = findddl(S,'all',nodeoutput,'free');

%% modele amortissement proportionnel
alpha = @(w) 0.01;
beta = @(w) 0.01;

%% affichage modele
figure(1);
clf
plot(S);

h1 = plotfacets(S,facesdirichlet,'EdgeColor','none','FaceColor','r','FaceAlpha',1);
h2 = plotfacets(S,facesefforts,'EdgeColor','none','FaceColor','y','FaceAlpha',1);
h3 = plotridges(S,numridges,'MarkerSize',20);
legend([h1,h2,h3],'Dirichlet','Neumann','Output')

%% calcul vibrations
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
Vz = V(findddl(S,'UZ'),:);

figure(10)
clf
scanmode = 1:6;
nbm = length(scanmode); 
for k=scanmode
    subplot(ceil(sqrt(nbm)),ceil(nbm/ceil(sqrt(nbm))),find(k==scanmode))
    title(['w_' num2str(k) ' = '  num2str(wres(k))]);
    scal = 1/max(max(abs(V(:,k))))*max(max(abs(getcoord(getnode(S)))))*0.4;
    plot(FENODEFIELD(sqrt(Vx(:,k).^2+Vy(:,k).^2+Vz(:,k).^2)),S+scal*V(:,k));
    pause(0.5)
end
