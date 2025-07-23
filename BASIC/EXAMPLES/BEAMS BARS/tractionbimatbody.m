a1=1;
a2=2;
a=1/2;
n=2*1;
h=1/n;
mat=MATERIALS();
mat{1}=ELAS_ISOT('E',a1,'S',1);
mat{2}=ELAS_ISOT('E',a2,'S',1);

P = POINT([0;a;1]);
L1 = mesh(LINE(P(1),P(2)),n/2,mat{1});
L2 = mesh(LINE(P(2),P(3)),n/2,mat{2});

S = union(L1,L2) ;
S=final(S,'norenum');
S=addcl(S,P(1),'UX',0);
K = calc_rigi(S);
f=bodyload(S,S,'FX',1);

u=K\f;

uh = eval_sol(S,u,P(3),'UX');


% u1 = inline('(x-x.^2/2)/a2+(1/a1-1/a2)*(a-a.^2/2)','a','a1','a2','x');
u1 = @(a,a1,a2,x) (x-x.^2/2)/a2+(1/a1-1/a2)*(a-a.^2/2);

uex = u1(a,a1,a2,1)
fprintf('Erreur en deplacement = %3d\n',abs(uh-uex)/abs(uex))

