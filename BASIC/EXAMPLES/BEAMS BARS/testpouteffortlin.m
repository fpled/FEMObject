r=10;
L=1;
d=0.1;
E=70e9;
I=pi*d^4/64;
Sec = pi*d^2/4;
mat = ELAS_BEAM('E',E,'NU',0.3,'RHO',7600,'S',Sec,'IZ',I,'IY',I,'IX',I*2);

P1=POINT([0,0]);
P2=POINT([L,0]);
P3=POINT([L/2,0]);
S1=mesh(LIGNE(P1,P2),r);
S=MODEL('PLAN');
S=addelem(S,'BEAM',S1,'mat',mat);
S=final(S,'norenum');

S=addcl(S,P1,{'U','R'},0);

K=calc_rigi(S);

ampl=2e6; 
% fun = inline('ones(size(x,1),1)','x');
fun = @(x) ones(size(x,1),1);
f=bodyload(S,LIGNE(P1,P2),{'FY'},fun);

figure(1)
q=K\f;
q=unfreevector(S,q);
plot(S+ampl*q)
x=getcoord(S.node);
x=x(:,1);

a=1/2*(x.^4/12-x.^3*L/3+L^2*x.^2/2)/E/I;
b=q(2:3:end);
fprintf('effort vertical constant\n')
norm(a-b)/norm(a)

% fun = inline('[x(:,1),ones(size(x,1),1)]','x');
fun = @(x) [x(:,1),ones(size(x,1),1)];
f=bodyload(S,LIGNE(P1,P2),{'FX'},fun);
ampl=1e9; 
figure(1)
q=K\f;
q=unfreevector(S,q);
plot(S+ampl*q)
a=1/2*(L^2*x-x.^3/3)/E/Sec;
b=q(1:3:end,1);
fprintf('effort longitudinal affine\n')
norm(a-b)/norm(a)

% fun = inline('x(:,1)','x');
fun = @(x) x(:,1);
f=bodyload(S,LIGNE(P1,P2),{'MZ'},fun);
ampl=1e6; 
figure(1)
q=K\f;
q=unfreevector(S,q);
plot(S+ampl*q,'Color','r');

% fun = inline('x(:,1)','x');
fun = @(x) x(:,1);
f=nodalload(S,P2,{'FY'},1);
ampl=1e6; 
figure(1)
q=K\f;
q=unfreevector(S,q);
plot(S+ampl*q,'Color','b')
a=(L*x.^2/2-x.^3/6)/E/I;
b=q(2:3:end,1);
fprintf('effort vertical affine affine\n')
norm(a-b)/norm(a)

%s=calc_sigma(q,M);
