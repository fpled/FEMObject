

r1=20;r2=20;
tf = 1e2;
T = TIMEMODEL(0,tf,100);
k=46;

close all
display_=0;
L=1;H=1;
%k=46;
rho=7800;
C=460;
C=0.1;
mat=FOUR_ISOT('c',rho*C,'k',k);

P1=POINT([0,0]);
P2=POINT([L,0]);
P3=POINT([L,H]);
P4=POINT([0,H]);
P5=POINT([0,H/2]);
P6 = POINT([L/2,0]);
S=mesh(DOMAIN(2,P1,P3),r1,r2,mat);
%S=setmaterial(S,mat);
S=final(S);


D1  = STRAIGHTLINE(P1,P4);
D2  = STRAIGHTLINE(P2,P3);
D3  = STRAIGHTLINE(P1,P2);
D4  = STRAIGHTLINE(P4,P3);
L1  = LINE(P1,P5);
L2  = LINE(P1,P6);

%S=addcl(S,D2,'T',0);
S=addcl(S,D4,'T',0);
%S=addcl(S,D3,'T',0);


K=calc_rigi(S,mat);
M=calc_mass(S,mat);

% fun1 = inline('x(:,2)','x');
% fun2 = inline('-x(:,2)','x');
f1=surfload(S,L1,{'QN'},1e4);
f2=surfload(S,L2,{'QN'},1e4);
f=f1+f2;
q=K\f;

figure(10)
clf
plot(q,S)
title('solution stationnaire')


loadfun = @(N) rampe(N,0,tf);

%% EULERSOLVER explicit
N = EULERTIMESOLVER(T,'eulertype','explicit');
ft = f*loadfun(N);

qtexp = dsolve(N,ft,M,K);
figure(2)
clf
evol(qtexp,S);

%% EULERSOLVER implicit

N = EULERTIMESOLVER(T,'eulertype','implicit');
ft = f*loadfun(N);

qtimp = dsolve(N,ft,M,K);
figure(3)
clf
%qt = setevolparam(qt,'pausetime',0.2);
evol(qtimp,S);

%% DGTIMESOLVER (discontinuous galerkin in time)
N = DGTIMESOLVER(T,1); % DG degree 1
ft = f*loadfun(N);
qtdg = dsolve(N,ft,M,K);
figure(4)
clf
evol(qtdg,S);

