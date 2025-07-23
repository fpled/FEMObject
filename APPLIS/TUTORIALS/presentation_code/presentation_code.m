%% GEOMETRIE ET MAILLAGE
%--------------------------
D = DOMAIN(2,[0,0],[1,1]);
figure(1)
subplot(1,2,1)
axis off
axis image
title('geometrical object: DOMAIN')
plot(D)
S = mesh(D,6,6);
subplot(1,2,2)
plot(S)
title('mesh of DOMAIN')

%%

figure(2)
plot(S,'numnode','numelem')
title('display numbers of nodes and elements')

figure(3)
plot(S,'node')
title('display nodes of a mesh')

%% Manipulation maillages
P = POINT([0,0;1,2;1,1;2,2])
figure(4)
plot(P,'b*')
title('objet geometrique: POINT')

S1 = mesh(DOMAIN(2,P(1),P(2)),10,20);
S2 = mesh(DOMAIN(2,P(3),P(4)),10,10);
S = union(S1,S2);
%S = convertelem(S,'TRI3');
figure(5)
clf
plot(S,'group')
title('groupes d''elements')

%% Importation de maillage depuis castem
% Scast = cast2matlab_model('PLANE','C:\PROGS\CASTEM\EXAMPLES\Lshape_domain.txt');
% figure(70)
% plot(Scast);

%% Importation de maillage depuis gmsh
% can be a .geo file or a .msh file
Sgmsh=gmsh2femobject(2,[getfemobjectoptions('path') 'APPLIS/TUTORIALS/evolution_problems/example_canister.geo'],2);
figure(70)
plot(Sgmsh);

%%
B = create_boundary(S);
figure(6)
clf
plot(S)
plot(B,'node','Color','r')
title('boundary of a mesh (MODEL)')

L1 = LINE(P(4),POINT([2,1]));
L2 = STRAIGHTLINE(P(1),POINT([1,0]));
figure(7)
clf
plot(S)
plot(intersect(B,L1),'node','Color','g')
plot(intersect(B,L2),'node','Color','r')

title('Intersection of boundary and geometrical object')

%% MATERIALS
MAT = MATERIALS();
MAT{1} = ELAS_ISOT('E',1,'NU',0.3,'RHO',1);
MAT{2} = ELAS_ISOT('E',2,'NU',0.3,'RHO',1);

S = setmaterial(S,MAT{1},1);
S = setmaterial(S,MAT{2},2);

figure(8)
plotparamelem(S,'material')
title('Materiaux')

S = setmaterial(S,MAT{1});

%% FINALIZATION OF THE MODEL : construction of dof structure
S = final(S);

%% APPLY BOUNDARY CONDITIONS
S = addcl(S,L2,'U');
S = addcl(S,LINE(P(2),P(4)),'UY');
K = calc_rigi(S);

f = surfload(S,L1,'FX',1);

figure(100)
clf
plot(S)
hold on
vectorplot(S,'F',f)

%% STATIQUE
u = K\f ;

figure(101)
clf
subplot(1,2,1)
plot(S)
subplot(1,2,2)
plot(S+0.03*u,'Color','r')

s = calc_sigma(S,u);
figure(102)
clf
subplot(1,3,1)
plot(s,S+0.03*u,'compo','SMXX')
subplot(1,3,2)
plot(s,S+0.03*u,'compo','SMYY')
subplot(1,3,3)
plot(s,S+0.03*u,'compo','SMXY')


%% CALCUL 3D
test3D

%% CALCUL DE POUTRES
testportiquejc(5)
%%
testbigportique(3)

%% CALCUL DE COQUES
shell3D

%% THERMIQUE INSTATIONNAIRE
r=20;
testthermique_evol

%% DYNAMIQUE
test2Ddynamic

