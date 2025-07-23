%% Homogeneisation thermique 2D, materiaux aleatoires
dim=2;
size_img=[1 1];
r=50;
nelem=[r r];

D=DOMAIN(dim,[0 0],size_img);
V=getvolume(D);
S=mesh(D,nelem(1),nelem(2));
S=createddlnode(S,DDL('T'));

lst=LSCIRCLE(0.5,0.5,0.3338);
lste=double(lseval(lst,S));

kf=10;
km=1;
indicf=double(lste<0);

%% Introduction des variables aleatoires
rvm=RVUNIFORM(km-0.1,km+0.3);
rvf=RVNORMAL(kf,2);
R=RANDVARS(rvm,rvf);
[X,PC]=PCTPMODEL(R,'order',5,'pcg');

a1=setfree(BILINFORM(1,1,1),0);
A1=a1{S}(:,:);
a2=setfree(BILINFORM(1,1,indicf,0),0);
A2=a2{S}(:,:);
A=A1*X{1}+A2*(X{2}-X{1});

%% Conditions aux limites de Dirichlet
% Mode 1
S1e=S;
S1e = addcl(S1e,getedge(D,1),'T',@(x) x(:,1));
S1e = addcl(S1e,getedge(D,3),'T',@(x) x(:,1));
S1e = addcl(S1e,getedge(D,2),'T',@(x) x(:,1));
S1e = addcl(S1e,getedge(D,4),'T',@(x) x(:,1));

A1f=freematrix(A1,S1e);
A2f=freematrix(A2,S1e);

Af=A1f*X{1}+A2f*(X{2}-X{1});
Asep=SEPMATRIX(Af);

b1=-calc_nonhomogeneous_vector(S1e,A1);
b2=-calc_nonhomogeneous_vector(S1e,A2);
b=b1*X{1}+b2*(X{2}-X{1});
bsep=SEPVECTOR(b);

PGD = SEPSOLVER(getdim(Asep),'maxorder',30,'tol',1e-5,'display',true,...
    'maxiter',10,'updatedim',2:getdim(Asep),'update',1);
[u1sep,result1ref] = solve(Asep,bsep,PGD);
u1e=PCTPMATRIX(u1sep,PC,2:3);

% Mode 2
S2e=S;
S2e = addcl(S2e,getedge(D,1),'T',@(x) x(:,2));
S2e = addcl(S2e,getedge(D,3),'T',@(x) x(:,2));
S2e = addcl(S2e,getedge(D,2),'T',@(x) x(:,2));
S2e = addcl(S2e,getedge(D,4),'T',@(x) x(:,2));

b1=-calc_nonhomogeneous_vector(S2e,A1);
b2=-calc_nonhomogeneous_vector(S2e,A2);
b=b1*X{1}+b2*(X{2}-X{1});
bsep=SEPVECTOR(b);

PGD = SEPSOLVER(getdim(Asep),'maxorder',30,'tol',1e-5,'display',true,...
    'maxiter',10,'updatedim',2:getdim(Asep),'update',1);
[u2sep,result2ref] = solve(Asep,bsep,PGD);
u2e=PCTPMATRIX(u2sep,PC,2:3);

%% Conditions aux limites periodiques
S1p=S;
S1p=addcl(S1p,POINT([0 0; size_img(1) 0;0 size_img(2);size_img(1) size_img(2)]),'T');
S1p=addclperiodic(S1p,getedge(D,1),getedge(D,3),'T');
S1p=addclperiodic(S1p,getedge(D,2),getedge(D,4),'T');

A1f=freematrix(A1,S1p);
A2f=freematrix(A2,S1p);
Af=A1f*X{1}+A2f*(X{2}-X{1});
Asep=SEPMATRIX(Af);

x=getcoord(getnode(S)); % Contient les 2 macro modes
y=x(:,2);
x(:,2)=[];

% Mode 1
b1=-A1*x;
b1=freevector(S1p,b1);
b2=-A2*x;
b2=freevector(S1p,b2);
b=b1*X{1}+b2*(X{2}-X{1});
bsep=SEPVECTOR(b);

PGD = SEPSOLVER(getdim(Asep),'maxorder',30,'tol',1e-5,'display',true,...
    'maxiter',10,'updatedim',2:getdim(Asep),'update',1);
[u1sep,result3ref] = solve(Asep,bsep,PGD);
u1p=PCTPMATRIX(u1sep,PC,2:3);

% Mode 2
b1=-A1*y;
b1=freevector(S1p,b1);
b2=-A2*y;
b2=freevector(S1p,b2);
b=b1*X{1}+b2*(X{2}-X{1});
bsep=SEPVECTOR(b);

PGD = SEPSOLVER(getdim(Asep),'maxorder',30,'tol',1e-5,'display',true,...
    'maxiter',10,'updatedim',2:getdim(Asep),'update',1);
[u2sep,result4ref] = solve(Asep,bsep,PGD);
u2p=PCTPMATRIX(u2sep,PC,2:3);

%% Conditions aux limites de Neumann
S1n=mesh(D,nelem(1),nelem(2));
S1n=createddlnode(S1n,DDL('T'),DDL('QN'));
point=size_img/2;
S1n=addcl(S1n,POINT(point),'T',0);

A1f=freematrix(A1,S1n);
A2f=freematrix(A2,S1n);
Af=A1f*X{1}+A2f*(X{2}-X{1});
Asep=SEPMATRIX(Af);

P=cell(4,1);
L=cell(4,1);

P{1}=POINT([0 0]);
P{2}=POINT([size_img(1) 0]);
P{3}=POINT([size_img(1) size_img(2)]);
P{4}=POINT([0 size_img(2)]);

L{1}=LINE(P{1},P{2});
L{2}=LINE(P{2},P{3});
L{3}=LINE(P{3},P{4});
L{4}=LINE(P{4},P{1});

f01=calc_nonhomogeneous_vector(S1n,A1);
f02=calc_nonhomogeneous_vector(S1n,A2);
f0=f01*X{1}+f02*(X{2}-X{1});
f0=SEPVECTOR(f0);

% Mode 1
f1=surfload(S1n,L{2},'QN',1);
f1=f1+surfload(S1n,L{4},'QN',-1);
f1=f1*one(PC);
f1=SEPVECTOR(f1);

bsep=f1-f0;
PGD = SEPSOLVER(getdim(Asep),'maxorder',30,'tol',1e-5,'display',true,...
    'maxiter',10,'updatedim',2:getdim(Asep),'update',1);
[u1sep,result5ref] = solve(Asep,bsep,PGD);

% Mode 2
f2=surfload(S1n,L{1},'QN',-1);
f2=f2+surfload(S1n,L{3},'QN',1);
f2=f2*one(PC);
f2=SEPVECTOR(f2);

bsep=f2-f0;
PGD = SEPSOLVER(getdim(Asep),'maxorder',30,'tol',1e-5,'display',true,...
    'maxiter',10,'updatedim',2:getdim(Asep),'update',1);
[u2sep,result6ref] = solve(Asep,bsep,PGD);

u1n=PCTPMATRIX(u1sep,PC,2:3);
u2n=PCTPMATRIX(u2sep,PC,2:3);

%% Determination des conductivites homogeneisees
% Utilisation du taux d'entropie

nbreal=50;

kmat=zeros(nbreal,1);
kfib=zeros(nbreal,1);

Ke=zeros(nbreal,4);
Kp=zeros(nbreal,4);
Kn=zeros(nbreal,4);

kk=zeros(2);

for i=1:nbreal
    [a,r]=random(A);

    kmat(i)=randomeval(X{1},r);
    kfib(i)=randomeval(X{2},r);

    u1=randomeval(u1e,r);
    u1=unfreevector(S1e,u1);
    u2=randomeval(u2e,r);
    u2=unfreevector(S2e,u2);
    Ke(i,1)=u1'*a*u1;
    Ke(i,4)=u2'*a*u2;
    Ke(i,3)=u1'*a*u2;
    Ke(i,2)=Ke(i,3);

    u1=randomeval(u1p,r);
    u1=x+unfreevector(S1p,u1);
    u2=randomeval(u2p,r);
    u2=y+unfreevector(S1p,u2);
    Kp(i,1)=u1'*a*u1;
    Kp(i,4)=u2'*a*u2;
    Kp(i,3)=u1'*a*u2;
    Kp(i,2)=Kp(i,3);

    u1=randomeval(u1n,r);
    u1=unfreevector(S1n,u1);
    u2=randomeval(u2n,r);
    u2=unfreevector(S1n,u2);
    kk(1,1)=u1'*a*u1;
    kk(2,2)=u2'*a*u2;
    kk(1,2)=u1'*a*u2;
    kk(2,1)=kk(1,2);
    kk=inv(kk);
    Kn(i,:)=kk(:);
end

Ke=Ke/V;
Kp=Kp/V;
Kn=Kn/V;

K=zeros(nbreal,3);
K(:,1)=Ke(:,1);
K(:,2)=Kp(:,1);
K(:,3)=Kn(:,1);

disp([kmat kfib K])

