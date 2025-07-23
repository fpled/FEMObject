
L=1;H=1;

mat=ELAS_ISOT('E',1,'NU',0.3,'DIM3',1,'RHO',1);

P1=POINT([0,0]);
P2=POINT([L,0]);
P3=POINT([L,H]);
P4=POINT([0,H]);
r1=50;
r2=50;
S1=mesh(DOMAIN(2,P1,P3),r1,r2,'material',mat,'option','CONT');
S = S1 ;
S=convertelem(S,'TRI3'); %avec les TRI3 on n'a pas la solution exacte aux
%noeuds
S=final(S);


D1  = STRAIGHTLINE(P1,P4);
D2  = STRAIGHTLINE(P2,P3);
D3  = STRAIGHTLINE(P1,P2);

S=addcl(S,D1,'U',[0;0]);
S=addcl(S,D2,'U',[1;0]);

K=calc_rigi(S);
u0 = calc_init_dirichlet(S);
f = calc_fint(S,u0);
q=K\(-f);
q = unfreevector(S,q)+u0;

figure(1)
clf
plot(S,'Color','k')
plot(S+q,'Color','r')
s=calc_sigma(S,q);
figure(2)
clf
plot(s,S+q,'compo','SMXX')
colorbar
%plot(S,'FaceColor','w','EdgeColor','k')
%plot(S+q,'FaceAlpha',0.3,'FaceColor','r','EdgeColor','r')

