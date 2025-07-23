%% Homogeneisation thermique 3D
dim=3;
size_img=[1 1 1];
nelem=[30 30 3];

D=DOMAIN(dim,[0 0 0],size_img);
V=getvolume(D);
S=mesh(D,nelem(1),nelem(2),nelem(3));
S=createddlnode(S,DDL('T'));

lst=LSCYLINDER(3,0.5,0.5,0,0,0,1,0.3338);
lst=double(lseval(lst,S));
kf=10;
km=1;
% delta=1/50;
% indicf=1-(1+tanh(lst/delta))/2;
indicf=double(lst<0);
k=km+(kf-km)*indicf;

%% Construction de l'operateur
a=setfree(BILINFORM(1,1,k(:),0),0);
A=a{S}(:,:);

%% Conditions aux limites EBC
S1e=S;
S1e = addcl(S1e,getface(D,1),'T',@(x) x(:,1));
S1e = addcl(S1e,getface(D,3),'T',@(x) x(:,1));
S1e = addcl(S1e,getface(D,2),'T',@(x) x(:,1));
S1e = addcl(S1e,getface(D,4),'T',@(x) x(:,1));
S1e = addcl(S1e,getface(D,5),'T',@(x) x(:,1));
S1e = addcl(S1e,getface(D,6),'T',@(x) x(:,1));
Af=freematrix(A,S1e);
b=-calc_nonhomogeneous_vector(S1e,A);
u1ef=Af\b;
u1e=unfreevector(S1e,u1ef);
figure('Name','EBC - Mode 1')
plot_sol(S1e,u1ef,'EdgeColor','none')

S2e=S;
S2e = addcl(S2e,getface(D,1),'T',@(x) x(:,2));
S2e = addcl(S2e,getface(D,3),'T',@(x) x(:,2));
S2e = addcl(S2e,getface(D,2),'T',@(x) x(:,2));
S2e = addcl(S2e,getface(D,4),'T',@(x) x(:,2));
S2e = addcl(S2e,getface(D,5),'T',@(x) x(:,2));
S2e = addcl(S2e,getface(D,6),'T',@(x) x(:,2));
Af=freematrix(A,S2e);
b=-calc_nonhomogeneous_vector(S2e,A);
u2ef=Af\b;
u2e=unfreevector(S2e,u2ef);
figure('Name','EBC - Mode 2')
plot_sol(S2e,u2ef,'EdgeColor','none')

S3e=S;
S3e = addcl(S3e,getface(D,1),'T',@(x) x(:,3));
S3e = addcl(S3e,getface(D,3),'T',@(x) x(:,3));
S3e = addcl(S3e,getface(D,2),'T',@(x) x(:,3));
S3e = addcl(S3e,getface(D,4),'T',@(x) x(:,3));
S3e = addcl(S3e,getface(D,5),'T',@(x) x(:,3));
S3e = addcl(S3e,getface(D,6),'T',@(x) x(:,3));
Af=freematrix(A,S3e);
b=-calc_nonhomogeneous_vector(S3e,A);
u3ef=Af\b;
u3e=unfreevector(S3e,u3ef);
figure('Name','EBC - Mode 3')
plot_sol(S3e,u3ef,'EdgeColor','none')


% Determination du tenseur de conductivite homogeneise
Ke=zeros(3);
Ke(1,1)=u1e'*A*u1e;
Ke(1,2)=u1e'*A*u2e;
Ke(1,3)=u1e'*A*u3e;
Ke(2,1)=u2e'*A*u1e;
Ke(2,2)=u2e'*A*u2e;
Ke(2,3)=u2e'*A*u3e;
Ke(3,1)=u3e'*A*u1e;
Ke(3,2)=u3e'*A*u2e;
Ke(3,3)=u3e'*A*u3e;
Ke=Ke/V;

%% Conditions aux limites PBC
S1p=S;
c1=[0 size_img(1)];
c2=[0 size_img(2)];
c3=[0 size_img(3)];
for i=1:2
    for j=1:2
        for k=1:2
            S1p=addcl(S1p,POINT([c1(i) c2(j) c3(k)]),'T');
        end
    end
end
S1p = addclperiodic(S1p,getface(D,1),getface(D,3),'T');
S1p = addclperiodic(S1p,getface(D,2),getface(D,4),'T');
S1p = addclperiodic(S1p,getface(D,5),getface(D,6),'T');

x=getcoord(getnode(S)); % Contient les 3 macro modes
y=x(:,2);
z=x(:,3);
x(:,2:3)=[];

b=-A*x;
b=freevector(S1p,b);
Af=freematrix(A,S1p);
u1pf=Af\b;
u1p=unfreevector(S1p,u1pf);
u1pt=x+u1p;
figure('Name','PBC - Mode 1')
plot_sol(S,u1pt,'EdgeColor','none');

b=-A*y;
b=freevector(S1p,b);
u2pf=Af\b;
u2p=unfreevector(S1p,u2pf);
u2pt=y+u2p;
figure('Name','PBC - Mode 2')
plot_sol(S,u2pt,'EdgeColor','none');

b=-A*z;
b=freevector(S1p,b);
u3pf=Af\b;
u3p=unfreevector(S1p,u3pf);
u3pt=z+u3p;
figure('Name','PBC - Mode 3')
plot_sol(S,u3pt,'EdgeColor','none');

% Determination du tenseur de conductivite homogeneise
Kp=zeros(3);
Kp(1,1)=u1pt'*A*u1pt;
Kp(1,2)=u1pt'*A*u2pt;
Kp(1,3)=u1pt'*A*u3pt;
Kp(2,1)=u2pt'*A*u1pt;
Kp(2,2)=u2pt'*A*u2pt;
Kp(2,3)=u2pt'*A*u3pt;
Kp(3,1)=u3pt'*A*u1pt;
Kp(3,2)=u3pt'*A*u2pt;
Kp(3,3)=u3pt'*A*u3pt;
Kp=Kp/V;

%% Conditions aux limites NBC
S1n=mesh(D,nelem(1),nelem(2),nelem(3));
S1n=createddlnode(S1n,DDL('T'),DDL('QN'));
point=size_img/2;
S1n=addcl(S1n,POINT(point),'T',0);

P=cell(8,1);
H=cell(6,1);

x=size_img(1);
y=size_img(2);
z=size_img(3);

P{1}=POINT([0 0 0]);
P{2}=POINT([0 0 z]);
P{3}=POINT([0 y 0]);
P{4}=POINT([0 y z]);
P{5}=POINT([x 0 0]);
P{6}=POINT([x 0 z]);
P{7}=POINT([x y 0]);
P{8}=POINT([x y z]);

H{1}=PLANE(P{1},P{2},P{3});
H{4}=PLANE(P{5},P{6},P{7});

H{2}=PLANE(P{1},P{2},P{5});
H{5}=PLANE(P{3},P{4},P{7});

H{3}=PLANE(P{1},P{3},P{5});
H{6}=PLANE(P{2},P{4},P{6});

f=cell(3,1);

for i=1:dim
    f{i}=surfload(S1n,H{i},'QN',-1)+surfload(S1n,H{i+3},'QN',1);
end

f0=calc_nonhomogeneous_vector(S1n,A);

Af=freematrix(A,S1n);

b=f{1}-f0;
u1nf=solvesingular(Af,b);
u1n=unfreevector(S1n,u1nf);
figure('Name','NBC - Mode 1')
plot_sol(S1n,u1nf,'EdgeColor','none');

b=f{2}-f0;
u2nf=solvesingular(Af,b);
u2n=unfreevector(S1n,u2nf);
figure('Name','NBC - Mode 2')
plot_sol(S1n,u2nf,'EdgeColor','none');

b=f{3}-f0;
u3nf=solvesingular(Af,b);
u3n=unfreevector(S1n,u3nf);
figure('Name','NBC - Mode 3')
plot_sol(S1n,u3nf,'EdgeColor','none');

% Determination du tenseur de conductivite homogeneise
Kn=zeros(3);
Kn(1,1)=u1n'*A*u1n;
Kn(1,2)=u1n'*A*u2n;
Kn(1,3)=u1n'*A*u3n;
Kn(2,1)=u2n'*A*u1n;
Kn(2,2)=u2n'*A*u2n;
Kn(2,3)=u2n'*A*u3n;
Kn(3,1)=u3n'*A*u1n;
Kn(3,2)=u3n'*A*u2n;
Kn(3,3)=u3n'*A*u3n;
Kn=inv(Kn);


