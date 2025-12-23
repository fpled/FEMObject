% filename = 'trapezegros.geo';
filename = 'trapezegros.msh';
Sini = gmsh2femobject(3,fullfile(getfemobjectoptions('path'),'LEGACY','GSD','EXAMPLES','GROS',filename));
savepath = fileparts(mfilename('fullpath'));
%%
S = Sini;
facesefforts = [4];
% facesefforts = [26,31];
facesdirichlet = [5];

figure(1);
clf;
plotfacets(S,'EdgeColor','w','FaceColor','b','FaceAlpha',1,'legend',false)
light('Position',[3 0 3])
myprint(savepath,['trapeze3D'],'epsc2')
myprint(savepath,['trapeze3D'],'jpeg','300')

figure(2);clf;
plotridges(S,'EdgeColor','k','FaceColor','k','EdgeAlpha',1,'FaceAlpha',1,'legend',false)
h1=plotfacets(S,facesdirichlet,'EdgeColor','none','FaceColor','r','FaceAlpha',1)
h2=plotfacets(S,facesefforts,'EdgeColor','none','FaceColor','y','FaceAlpha',1)
legend([h1(1),h2(1)],'Dirichlet','Neumann')

%%
S = Sini;
mat = MATERIALS(ELAS_ISOT('E',1,'NU',0.3,'RHO',1));
S = setmaterial(S,mat{1});
S = splitgroupelem(S,3000);
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
K0 = calc_rigi(S);
E0 = 1;
% setfemobjectoptions('tolerancepoint',1e-7)

%% calcul deterministe
mat = MATERIALS(ELAS_ISOT('E',1,'NU',0.3,'RHO',1));
S = setmaterial(S,mat{1});

q = K0\f;
ampl = 0.01;
s = calc_sigma(S,q);

figure(3)
clf
plot_sol(S,q,'sigma',1,'ampl',1/500)

%%
m = 8;
lc = 0.5;

x = getcoord(getnode(S));

D = mesh(LINE(0,2),100);
[V,L] = KL(RFGAUSSIANEXP2(1,0.2,lc),m,D,'modes');
fprintf('\nmode #')
for i=1:m
    Emode{i} = L(i+1,i+1)*sin(i*pi/2*x(:,1));    
    mat = MATERIALS(ELAS_ISOT('E',FENODEFIELD(Emode{i}),'NU',0.3,'RHO',1));
    S = setmaterial(S,mat{1});    
    Kmode{i} = calc_rigi(S);
    fprintf('%d/',i);
end

%%
p = 6;
RV = RANDVARS(RVNORMAL(),m);
X = PCMODEL(RV,'order',1,'pcg');
PCX = getPC(X);
PC = POLYCHAOS(RANDPOLYS(PCX),p);
PCX = calc_masse(PCX,PC);


Ksto = K0; 
E = E0*ones(S.nbnode,1);
for i=1:m
    Ksto = Ksto + X{i}*Kmode{i};
    E = E + X{i}*Emode{i};
end

mat = MATERIALS(ELAS_ISOT('E',FENODEFIELD(E),'NU',0.3,'RHO',1));
S = setmaterial(S,mat{1});

clear Kmode
clear Emode

%% 
tic
fprintf('debut resolution pcg ... ')
% qpc = pcg(Ksto,f*one(PC),1e-4,200);
qpc = mypcg(Ksto,f*one(PC),'tol',1e-4,'noreortho','display');
fprintf('fin resolution pcg\n')
toc


%%
GSD = GSDSOLVER('tol',1e-4,'nbfoncmax',8,...
    'display',true,'errorindicator','none');
GSD = setparam(GSD,'type','power','update',true);
GSD = setparam(GSD,'type','arnoldi','restart',0,'orthocrit',1e-12);

[qrad,result]=solve(GSD,Ksto,f*one(PC));
qradsd = spectral_decomposition(qrad,'display');

%%
[qr,xi]=random(qrad);
Sr = randomeval(S,xi); 

figure(3)
clf
plot_sol(Sr,qr,'sigma',1,'ampl',1/500)
colorbar
