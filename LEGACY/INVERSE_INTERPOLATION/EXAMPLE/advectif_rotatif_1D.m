%% Modele deterministe
clear all

dens=30;


P1=POINT([0 0]);
P2=POINT([1 1]);
P3=POINT([.4 .4]);
P4=POINT([.6 .6]);
P5=POINT([.7 .7]);
P6=POINT([.9 .9]);
P7=POINT([1 0]);
P8=POINT([0 1]);
P9=POINT([0.1 0.1]);
D0=mesh(DOMAIN(2,P1,P2) , dens, dens);
h=1/dens;

S=D0;
S=createddlnode(S,DDL('u'));
S=final(S);

% Cls PERIODIQUE
S = addclperiodic(S,LINE(P1,P8),...
                    LINE(P7,P2),'u');
S = addclperiodic(S,LINE(POINT([h 0]),POINT([1-h 0])),...
                    LINE(POINT([h 1]),POINT([1-h 1])),'u');

plot(S)

%% Operateur
reac=1;

B=BILINFORM(1,1);
B=setfree(B,0);
K=freematrix(S,B{S}(:,:));
B=BILINFORM(0,0);
B=setfree(B,0);
K=K+reac*freematrix(S,B{S}(:,:));

B  = BILINFORM(0,1,{1,0},0);
B  = setfree(B,0);
AX = freematrix(S,B{ S }(:,:));
B  = BILINFORM(0,1,{0,1},0);
B  = setfree(B,0);
AY = freematrix(S,B{ S }(:,:));

A_full={K,AX,AY};

% rhs
aa=0.05;
bb=0.1;
X=getcoord(getnode(S));
f= (bb-sqrt((sum((X-[ones(size(X,1),1)*0.25 , ones(size(X,1),1)*0.5] ).^2,2))) )/(bb-aa);
f(f<0)=0;
f(f>=1)=1;

aa=0.05;
bb=0.1;
X=getcoord(getnode(S));
f2= (bb-sqrt((sum((X-[ones(size(X,1),1)*0.75 , ones(size(X,1),1)*0.5] ).^2,2))) )/(bb-aa);
f2(f2<0)=0;
f2(f2>=1)=1;
f=f-f2;

plot(f,S)

L=LINFORM(0,f,0);

f=L{S}(:);
f=full(f);


%%

%
D=50;
teta=rand()*pi*2;

alpha= D*cos(teta);
beta = D*sin(teta);

A_xi= K + alpha*AX + beta*AY;

u=A_xi\f;
clf
plot(unfreevector(S,u),S)
colorbar

%% Finalisation

perm = symamd( A_xi );

A = cellfun(@(x) x(perm,perm),A_full,'UniformOutput',0);
b = full(f(perm));

PhiA  = @(x) [ ones(1,size(x,1)) ;  D*cos(pi*2*x');  D*sin(pi*2*x')] ;
dim_p = 1;



















