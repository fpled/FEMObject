r=20;
m=[1,2,3,4,5]
L=1;
d=0.1;

MAT = ELAS_BEAM('E',70e9,'NU',0.3,'RHO',7600,'S',pi*d^2/4,'IZ',pi*d^4/64,'IY',pi*d^4/64,'IX',pi*d^4/32);

P1=POINT([0,0]);
P2=POINT([L,0]);
P3=POINT([L/2,0]);
S1=mesh(LINE(P1,P2),r,MAT);
S=convertelem(S1,'BEAM','param',VECTOR([0;1]));
S=final(S);
S=addcl(S,P1,{'U','R'},0);
%S=addcl(S,P2,{'U','R'},0);

K=calc_rigi(S);
M=calc_mass(S);


ampl=1; 


[V,D]=calc_mode(K,M,1:max(m));
n=7; ampli = cos(linspace(0,pi,n))*ampl;
figure(1);
axis0 =[-L/5,L+L/5,-2*L,2*L];
varargout{1}=V;
varargout{2}=D;

for i=1:length(ampli) 
    clf;
    for k=m(:)'
        subplot(1,length(m),find(k==m))

        plot(S+ampli(i)*V(:,k));

        axis(axis0);
    end
%    FILM(i)=getframe;
    pause(0.01);
end
% movie(FILM,1,24)
ampl = 10000; 
f=nodalload(S,P2,'FY',100);    
q=K\f;
figure(2)
clf

plot(S+ampl*q)

figure(3)
clf
s=calc_sigma(S,q,'smooth');
plot(s(2),S);title('MOMZ')
