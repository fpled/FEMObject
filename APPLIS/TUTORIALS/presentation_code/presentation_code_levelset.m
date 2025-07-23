%% LEVELSET

S0 = mesh(DOMAIN(2,[0,0],[1,1]),15,15);
S0 = convertelem(S0,'TRI3')
ls = LSCIRCLE(.3,.4,.2);

figure(105)
plot(S0)
%plot(ls,S,'surface')
contourplot(ls,S0)

%% Manipulation levels-sets
figure(106)
clf
subplot(1,3,1)
ls1 = LSCIRCLE(.3,.3,.2)
plot(S0)
plot(ls1,S0)
contourplot(ls1,S0)
title('ls1')
subplot(1,3,2)
ls2 = LSELLIPSE(.4,.6,.2,.3,1,0)
plot(S0)
plot(ls2,S0)
contourplot(ls2,S0)
title('ls2')
subplot(1,3,3);
ls3 = complement(union(ls1,ls2));
plot(S0)
plot(ls3,S0)
contourplot(ls3,S0)
title('union de ls1 et ls2')

%% operations boolennes
% complement, union, interset, setdiff
%%

S0 = mesh(DOMAIN(2,[0,0],[1,1]),20,20);
omega = RVUNIFORM(1,4);
phi = RVUNIFORM(0,2*pi);
funa = @(x,a) cos(pi*x*a(1)+a(2))*0.1;
ls = LSROUGHLINE(0,0.8,0,1,funa,RANDVARS(omega,phi));
tic

%
S = LSMODEL(S0,random(ls));
toc
%S = convertelem(S,'TRI3');
tic
figure(1);clf;plot(getlevelset(S,1),S);
contourplot(getlevelset(S,1),S,'Color','k');plot(S)
figure(2)
clf
toc
tic
Sdiv = lsdivideelem(S);
toc
plot(Sdiv);


%% 
S0 = mesh(DOMAIN(2,[0,0],[1,1]),25,25);
S0 = convertelem(S0,'TRI3')

MAT = MATERIALS();
MAT{1} = ELAS_ISOT('E',1,'NU',0.3,'RHO',1);

S = LSMODEL(S0,ls3);
S = setmaterial(S,MAT{1});
S = final(S);
figure(1000)
clf
lsplot(S);

S = addcl(S,LINE([0,0],[0,1]),'U');
f = surfload(S,LINE([1,0],[1,1]),'FX',1);
K = calc_rigi(S);
u = K\f;
s = calc_sigma(S,u);

figure(1001)
clf
plot(s,S,'compo','SMXX')
contourplot(ls3,S)
figure(1002)
clf
plot_sol(S,.05*u,'selgroup',1:2)

%%
Sf = lsdivideelem(S)
Sf = final(Sf);
figure(1003)
clf
plot(Sf,'Color','r')
plot(S)
%uf = transfer(S,Sf,unfreevector(S,u))

%% VARIABLES ALEATOIRES

r = RVNORMAL(1,0.2);
pdfplot(r,'r')

random(r)

%%
rl = RVLOGNORMAL(0.2,0.2,0,'stat');

rlpc = decomppc(rl,'order',3,'nbgauss',10);
getPC(rlpc)
figure(1000)

pdfplot(rl)
hold on
pdfplot(rlpc,'r','nbs',1e5)

%% LEVELSET ALEATOIRE

RV = RANDVARS(RVUNIFORM(.2,.4),RVUNIFORM(.4,.7),RVUNIFORM(.1,.3));
S = mesh(DOMAIN(2,[0,0],[1,1]),15,15);
ls = LSCIRCLE(RV{1},RV{2},setnumber(rlpc,3))

figure
plot(S)
contourplot(random(ls),S)
%%
figure(1400)
clf
lsr = random(ls);
plot(S)
plot(lsr,S)
contourplot(lsr,S)







