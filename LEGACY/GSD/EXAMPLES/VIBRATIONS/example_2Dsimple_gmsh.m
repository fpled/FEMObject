%% importation du modele
Dom = DOMAIN(2,[0,0],[2,1]);
S = gmsh(Dom,0.1,'all');

%% construction du modele EF et CL
facesdirichlet = 4; 
facesefforts = 2; 

mat = MATERIALS(ELAS_ISOT('E',1,'NU',0.3,'RHO',1));
S = setmaterial(S,mat{1});
S = final(S);
S = addcl(S,getedge(Dom,facesdirichlet),'U');

f = surfload(S,getedge(Dom,facesefforts),'FY',1);

K1 = calc_rigi(S);
M1 = calc_mass(S);

%% output
numridges = 3;
nodeoutput = [];
for l=numridges
    nodeoutput = union(nodeoutput,getnumber(getnode(getridge(S,l))));
end
ddloutput = findddl(S,'all',nodeoutput,'free');

%% modele amortissement proportionnel
alpha = @(w) .01;
beta = @(w) 0.01;

%% Affichage modele
figure(1)
clf
plot(S);

h1 = plotfacets(S,facesdirichlet,'FaceColor','r','EdgeColor','none','FaceAlpha',1);
h2 = plotfacets(S,facesefforts,'FaceColor','y','EdgeColor','none','FaceAlpha',1);
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
