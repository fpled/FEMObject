%% importation du modele
% filename = 'trapezegros.geo';
filename = 'trapezegros.msh';
Sini = gmsh2femobject(3,fullfile(getfemobjectoptions('path'),'LEGACY','GSD','EXAMPLES','GROS',filename));
savepath = fileparts(mfilename('fullpath'));

%% construction du modele EF et CL
facesefforts = [3];
facesdirichlet = [5];

mat = MATERIALS(ELAS_ISOT('E',1,'NU',0.3,'RHO',1));
S = setmaterial(S,mat{1});
S = splitgroupelem(S,8000);
S = final(S,'renum');

nodesdirichlet=[];
for l=facesdirichlet % 20:23
    nodesdirichlet = union(nodesdirichlet,getnumber(getnode(getfacet(S,l))));
end
nodesefforts=[];
for l=facesefforts % 40
    nodesefforts = union(nodesefforts,getnumber(getnode(getfacet(S,l))));
end

S = addcl(S,nodesdirichlet,'U');    
f = nodalload(S,nodesefforts,'FY',-1);  

K1 = calc_rigi(S);
M1 = calc_mass(S);

% setfemobjectoptions('tolerancepoint',1e-7)

%% output
numpeaks = [3,15];
nodeoutput = [];
for l=numpeaks
    nodeoutput = union(nodeoutput,getnumber(getnode(getpeak(S,l))));
end
ddloutput = findddl(S,'all',nodeoutput,'free');

%% affichage modele
figure(1)
clf
plotfacets(S,'EdgeColor','w','FaceColor','b','FaceAlpha',1,'legend',false)
light('Position',[3  0  3])
myprint(savepath,['trapeze3D'],'epsc2')
myprint(savepath,['trapeze3D'],'jpeg','300')

figure(2)
clf
plotridges(S,'EdgeColor','k','FaceColor','k','EdgeAlpha',1,'FaceAlpha',1,'legend',false);
h1 = plotfacets(S,facesdirichlet,'EdgeColor','none','FaceColor','r','FaceAlpha',1);
h2 = plotfacets(S,facesefforts,'EdgeColor','none','FaceColor','y','FaceAlpha',1);
legend([h1(1),h2(1)],'Dirichlet','Neumann')

figure(3)
clf
plotridges(S,'EdgeColor','k','FaceColor','k','EdgeAlpha',1,'FaceAlpha',1,'legend',false)
h1 = plotfacets(S,facesdirichlet,'EdgeColor','none','FaceColor','r','FaceAlpha',1);
h2 = plotfacets(S,facesefforts,'EdgeColor','none','FaceColor','y','FaceAlpha',1);
h3 = plotpeaks(S,numpeaks,'MarkerSize',20);
legend([h1(1),h2(1),h3(1)],'Dirichlet','Neumann','Output')

%% modele amortissement proportionnel
alpha = @(w) 0.001;
beta = @(w) 0.001;

%% calcul vibrations sans amort
scanm = 1:30
[V,s] = calc_mode(K1,M1,scanm);
s = diag(s);
wres = sqrt(s)

mm = 8;
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
    plot(FENODEFIELD(sqrt(Vy(:,k).^2+Vx(:,k).^2)),S+scal*V(:,k));
    pause(0.5)
end
