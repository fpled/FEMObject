function varargout = test2Deffortsurf(r1,r2,varargin)
close all
L=3;H=1;

mat=ELAS_ISOT('E',1,'NU',0.3,'DIM3',1,'RHO',1);

P1=POINT([0,-H/2]);
P2=POINT([L,-H/2]);
P3=POINT([L,H/2]);
P4=POINT([0,H/2]);

S1=mesh(DOMAIN(2,P1,P3),r1,r2,'mat',mat);
S = S1 ;
S=setoption(S,'CONT');
S=final(S);


D1  = STRAIGHTLINE(P1,P4);
D2  = STRAIGHTLINE(P2,P3);
D3  = STRAIGHTLINE(P1,P2);
S=addcl(S,P1,'U',[0;0]);
S=addcl(S,P2,'UY',0);

K=calc_rigi(S);
% fun1 = inline('x(1,2)','x');
% fun2 = inline('-x(1,2)','x');
fun1 = @(x) x(1,2);
fun2 = @(x) -x(1,2);
plot(D1)
plot(D2)
f1=surfload(S,D1,{'FX'},fun1);
f2=surfload(S,D2,{'FX'},fun2);
f=f1+f2;

q=K\f;

q=unfreevector(S,q);

figure(2)
clf
varargout{1}=q;
varargout{2}=f;
[rep,pos]=ischarin('ampl',varargin);

if rep
    ampl = varargin{pos+1};
else
    ampl=0.1; 
end
plot(S+ampl*q)
rep=findddl(S,{'UX'},POINT([L,H/2]));
q(rep)
class(ampl*q)
s=calc_elemfield(S,'sigmanode',q);
figure(3)
plot(s,S+ampl*q,'compo','SMXX');
colorbar
