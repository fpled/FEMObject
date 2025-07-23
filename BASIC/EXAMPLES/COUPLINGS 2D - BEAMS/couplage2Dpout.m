warning('c''est un peu rafistole : le raccord entre modele 2D et poutre est meilleur avec des coques')


r1=20;
r2=5;
H=1;
L=3;
d=0.1;
mat=MATERIALS();
mat{1} = ELAS_ISOT('E',70e9,'NU',0.3,'DIM3',d,'RHO',7600);
mat{2} = ELAS_BEAM('E',70e9,'NU',0.3,'RHO',7600,'S',pi*d^2/4,'IZ',pi*d^4/64,'IY',pi*d^4/32,'IX',pi*d^4/32);

P1=POINT([0,0]);
P2=POINT([L,0]);
P3=POINT([L,H]);
P4=POINT([L,2*H]);
P5=POINT([0,2*H]);
P6=POINT([2*L,H]);
D1 = STRAIGHTLINE(P1,P5);

S1 = mesh(DOMAIN(2,P1,P4),r1,2*r2,mat{1});
S2 = mesh(LINE(P3,P6),r1,mat{2});
S2 = convertelem(S2,'BEAM');
S = union(S1,S2);

S=final(S);

S=addcl(S,D1,{'U'},[0;0]);
S=addcl(S,P3,'R',0);

K=calc_rigi(S);


F=nodalload(S,P6,'FY',1e5);

q=K\F;


figure(2)
clf
ampl=1;
plot(S,'Color','b')
plot(S+ampl*q)


figure(3)
clf
ampl=1;
s=calc_sigma(S,q);
plot(s,S+ampl*q);


%% transfert de la solution sur le modele 2D et affichage de la deformee amplifiee
q = unfreevector(S,q);
S1 = final(S1);
q1 = transfer(S,S1,q);

figure(10)
clf
ampl=1000;
s1 = calc_sigma(S1,q1);
subplot(1,3,1)
plot(s1,S1+ampl*q1,'COMPO','SMXX')
subplot(1,3,2)
plot(s1,S1+ampl*q1,'COMPO','SMXY')
subplot(1,3,3)
plot(s1,S1+ampl*q1,'COMPO','SMYY')
