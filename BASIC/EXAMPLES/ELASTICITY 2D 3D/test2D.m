r=5;
m=3;

e=0.2;L=1;H=1;
r1=r;r2=r*L/e;
mat=MATERIALS();
mat{1}=ELAS_ISOT('E',1,'NU',0.3,'DIM3',1,'RHO',1);
mat{2}=ELAS_ISOT('E',1,'NU',0.3,'DIM3',1,'RHO',1);
mat{3}=ELAS_ISOT('E',1,'NU',0.3,'DIM3',1,'RHO',1);

P1=POINT([0,0]);
P2=POINT([e,H]);
P3=POINT([e,H-e]);
P4=POINT([L-e,H]);
P5=POINT([L-e,0]);
P6=POINT([L,H]);
P7=POINT([L,0]);

S1=mesh(DOMAIN(2,P1,P2),r1,r2);
S2=mesh(DOMAIN(2,P3,P4),floor(r2*(L-2*e)/L),r1);
S3=mesh(DOMAIN(2,P5,P6),r1,r2);
S=union(S1,S2,S3);
%S=convertelem(S,'TRI3');
S=setmaterial(S,mat{1},1);
S=setmaterial(S,mat{1},3);
S=setmaterial(S,mat{1},2);
%S = concatgroupelem(S);
%S = splitgroupelem(S,2000);

%S=MODEL('PLANE');
%S=addelem(S,'QUA4',S1,'mat',mat{1});
%S=addelem(S,'QUA4',S2,'mat',mat{2});
%S=addelem(S,'QUA4',S3,'mat',mat{3});
%S=convertelem(S,'TRI3');

S=final(S);

D1  = STRAIGHTLINE(P1,P7);
D3  = STRAIGHTLINE(P6,VECTOR([1;0]));

S=addcl(S,D1,'U');

K=calc_rigi(S);


M=calc_mass(S);
[V,D]=calc_mode(K,M,1:m);
n=100; ampl = cos(linspace(0,2*pi,n))/10;
figure(1);
axis0 =[-L/3,L+L/3,-H/3,H+H/3];
for i=1:length(ampl) ;clf;plot(S+ampl(i)*V(:,m));axis(axis0);pause(0.01);end


f=surfload(S,D3,'FY',-1);

q=K\f;


figure(2)
clf
ampl=0.05;
plot(S+ampl*q)
pause(0.3)

s=calc_sigma(S,q,'node');
sm=calc_sigma(S,q,'smooth');
e=calc_epsilon(S,q,'node');

figure(3)
subplot(2,2,1)
plot(s,S+ampl*q,'compo','SMXX')
cax0 = caxis;
subplot(2,2,2)
plot(sm,S+ampl*q,'compo','SMXX')
caxis(cax0)

subplot(2,2,3)
plot(s,S+ampl*q,'compo','SMYY')
cax0 = caxis;
subplot(2,2,4)
plot(sm,S+ampl*q,'compo','SMYY')
caxis(cax0)


