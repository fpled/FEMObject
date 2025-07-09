
L = 1;
H = 1;

mat = ELAS_ISOT('E',1,'NU',0,'DIM3',1,'RHO',1);

P1 = POINT([0,0]);
P2 = POINT([L,0]);
P3 = POINT([L,H]);
P4 = POINT([0,H]);
r1 = 60;
r2 = 60;

S1 = mesh(DOMAIN(2,P1,P3),r1,r2,'material',mat,'option','CONT');
S = S1 ;
S = convertelem(S,'TRI3'); %avec les TRI3 on n'a pas la solution exacte aux noeuds
S = final(S);

D1 = DROITE(P1,P4);
D2 = DROITE(P2,P3);
D3 = DROITE(P1,P2);

S = addcl(S,D1,'U',[0;0]);

K = calc_rigi(S);

f = bodyload(S,[],'FX',1);

q = K\f;

% u1 = inline('(L*x-x.^2/2)','L','x');
u1 = @(L,x) (L*x-x.^2/2);
uex = u1(L,L)
uapp = eval_sol(S,q,POINT([L,1]),'UX')
erroru = norm(uapp-uex)/norm(uex)

ampl = L/max(abs(q));

figure(1)
clf
plot(S,'facecolor','b','edgecolor','b','facealpha',0.3)
plot(S+ampl*q,'facecolor','r','edgecolor','r','facealpha',0.3)
