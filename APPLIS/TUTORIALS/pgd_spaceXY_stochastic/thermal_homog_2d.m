%% Homogeneisation thermique 2D
dim=2;
size_img=[1 1];
r=50;
nelem=[r r+2];

D=DOMAIN(dim,[0 0],size_img);
V=getvolume(D);
S=mesh(D,nelem(1),nelem(2));
S=createddlnode(S,DDL('T'));

lst=LSCIRCLE(0.5,0.5,0.3338);
lste=double(lseval(lst,S));
kf=10;
km=1;
% delta=1/50;
% indicf=1-(1+tanh(lst/delta))/2;
indicf=double(lste<0);
k=km+(kf-km)*indicf;
figure
surface(reshape(k,nelem(2)+1,nelem(1)+1))

%% Construction de l'operateur
a=setfree(BILINFORM(1,1,k(:),0),0);
A=a{S}(:,:);

%% Conditions aux limites de Dirichlet
% Mode 1
S1e=S;
S1e = addcl(S1e,getedge(D,1),'T',@(x) x(:,1));
S1e = addcl(S1e,getedge(D,3),'T',@(x) x(:,1));
S1e = addcl(S1e,getedge(D,2),'T',@(x) x(:,1));
S1e = addcl(S1e,getedge(D,4),'T',@(x) x(:,1));
Af=freematrix(A,S1e);
b=-calc_nonhomogeneous_vector(S1e,A);
u1ef=Af\b;
u1e=unfreevector(S1e,u1ef);
figure('Name','EBC - Mode 1')
plot_sol(S1e,u1ef,'EdgeColor','none')

% Mode 2
S2e=S;
S2e = addcl(S2e,getedge(D,1),'T',@(x) x(:,2));
S2e = addcl(S2e,getedge(D,3),'T',@(x) x(:,2));
S2e = addcl(S2e,getedge(D,2),'T',@(x) x(:,2));
S2e = addcl(S2e,getedge(D,4),'T',@(x) x(:,2));
Af=freematrix(A,S2e);
b=-calc_nonhomogeneous_vector(S2e,A);
u2ef=Af\b;
u2e=unfreevector(S2e,u2ef);
figure('Name','EBC - Mode 2')
plot_sol(S2e,u2ef,'EdgeColor','none')

%% Conditions aux limites periodiques
S1p=S;
S1p=addcl(S1p,POINT([0 0; size_img(1) 0;0 size_img(2);size_img(1) size_img(2)]),'T');
S1p=addclperiodic(S1p,getedge(D,1),getedge(D,3),'T');
S1p=addclperiodic(S1p,getedge(D,2),getedge(D,4),'T');

x=getcoord(getnode(S)); % Contient les 2 macro modes
y=x(:,2);
x(:,2)=[];

% Mode 1
b=-A*x;
b=freevector(S1p,b);
Af=freematrix(A,S1p);
u1pf=Af\b;
u1p=unfreevector(S1p,u1pf);
u1pt=x+u1p;
figure('Name','PBC - Mode 1')
plot_sol(S,u1pt,'EdgeColor','none');

% Mode 2
b=-A*y;
b=freevector(S1p,b);
u2pf=Af\b;
u2p=unfreevector(S1p,u2pf);
u2pt=y+u2p;
figure('Name','PBC - Mode 2')
plot_sol(S,u2pt,'EdgeColor','none');

%% Conditions aux limites de Neumann
S1n=mesh(D,nelem(1),nelem(2));
S1n=createddlnode(S1n,DDL('T'),DDL('QN'));
point=size_img/2;
S1n=addcl(S1n,POINT(point),'T',0);

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

f=calc_nonhomogeneous_vector(S1n,A);

Af=freematrix(A,S1n);

% Mode 1
f1=surfload(S1n,L{2},'QN',1);
f1=f1+surfload(S1n,L{4},'QN',-1);
b=f1-f;
u1nf=solvesingular(Af,b);
u1n=unfreevector(S1n,u1nf);
figure('Name','NBC - Mode 1')
plot_sol(S1n,u1nf,'EdgeColor','none');

% Mode 2
f2=surfload(S1n,L{1},'QN',-1);
f2=f2+surfload(S1n,L{3},'QN',1);
b=f2-f;
u2nf=solvesingular(Af,b);
u2n=unfreevector(S1n,u2nf);
figure('Name','NBC - Mode 2')
plot_sol(S1n,u2nf,'EdgeColor','none');

%% Determination des conductivites homogeneisees
% Utilisation du taux d'entropie
Ke=zeros(2);
Ke(1)=u1e'*A*u1e;
Ke(4)=u2e'*A*u2e;
Ke(2)=u1e'*A*u2e;
Ke(3)=Ke(2);
Ke=Ke/V;

Kp=zeros(2);
Kp(1)=u1pt'*A*u1pt;
Kp(4)=u2pt'*A*u2pt;
Kp(2)=u1pt'*A*u2pt;
Kp(3)=Kp(2);
Kp=Kp/V;

Kn=zeros(2);
Kn(1)=u1n'*A*u1n;
Kn(4)=u2n'*A*u2n;
Kn(2)=u1n'*A*u2n;
Kn(3)=Kn(2);
Kn=Kn/V;
Kn=inv(Kn);
