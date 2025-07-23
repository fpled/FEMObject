n=10;
P=POINT([0,0;1,0;1,1;0,1;1,1/4;1,3/4]);

MAT = ELAS_ISOT('E',1,'NU',0.,'RHO',1,'DIM3',1);
D1=LINE(P(1),P(2));
D2=LINE(P(2),P(3));
D3=LINE(P(3),P(4));
D4=LINE(P(4),P(1));
D5=LINE(P(5),P(6));

S=mesh(DOMAIN(2,P(1),P(3)),2*n,2*n,MAT);
S=final(S);
S=addcl(S,D4,'U',[0;0]);
K=calc_rigi(S);
M=calc_mass(S);

ampl=0.1;
f=surfload(S,D5,{'FX','FY'},[-1;0]);


nt=100;
T = 3;
t = linspace(0,T,nt+1);
dt=T/nt ;
alpha=0.05;
gamma = 1/2 +alpha ;
beta = (1+alpha)^2/4 ;
tc = T/4 ;
rep = find(t<=tc);repo = find(t>tc);
ft=f*[t(rep)/tc,ones(1,length(repo))];
ut = zeros(size(K,1),nt+1);
vt = zeros(size(K,1),nt+1);
at = zeros(size(K,1),nt+1);

disp('resolution temporelle')
R = chol(M+(1-alpha)*K*dt^2*beta);
at(:,1)=M\ft(:,1);
for i=1:nt
    pourcentage(i,nt)
    at(:,i+1) = R\(R'\((1-alpha)*ft(:,i+1)+alpha*ft(:,i)-alpha*K*ut(:,i)-(1-alpha)*K*(ut(:,i)+vt(:,i)*dt+(1/2-beta)*at(:,i)*dt^2)) );   
    vt(:,i+1) = vt(:,i) + dt*((1-gamma)*at(:,i)+gamma*at(:,i+1));
    ut(:,i+1) = ut(:,i) + dt*vt(:,i) + dt^2*((1/2-beta)*at(:,i)+beta*at(:,i+1));
end
fprintf('\n')
clear s

s=calc_sigma(S,ut,'smooth')  ;  
sm=s;
pas=1;
figure(2)
plot(s(:,end),S,'compo','SMXX');
cax=caxis;

figure(2)
for i=1:pas:nt
    clf
    plot(s(:,i),S,'compo','SMXX');
    title(['time ' num2str(t(i)) ' s'])
    caxis(cax)
    pause(0.01)
end

figure(3)
plot(sm(:,end),S,'compo','SMXX','surface');
ax=axis;
cax=caxis;
ax(1)=ax(1)-(ax(2)-ax(1))/2;
ax(2)=ax(2)+(ax(2)-ax(1))/2;
ax(3)=ax(3)-(ax(4)-ax(3))/2;
ax(4)=ax(4)+(ax(4)-ax(3))/2;

figure(3)
for i=1:pas:nt
    clf
    plot(sm(:,i),S,'compo','SMXX','surface');
    title(['time ' num2str(t(i)) ' s'])
    axis(ax)
    caxis(cax)
    axis on
    pause(0.01)
end

figure(4)
ut=unfreevector(S,ut);
rep = findddl(S,'UX');

plot(FENODEFIELD(ut(rep,end)),S,'surface');
ax=axis;
cax=caxis;
ax(1)=ax(1)-(ax(2)-ax(1))/2;
ax(2)=ax(2)+(ax(2)-ax(1))/2;
ax(3)=ax(3)-(ax(4)-ax(3))/2;
ax(4)=ax(4)+(ax(4)-ax(3))/2;
for i=1:pas:nt
    clf
    plot(FENODEFIELD(ut(rep,i)),S,'compo','SMXX','surface');
    title(['time ' num2str(t(i)) ' s'])
    axis(ax)
    caxis(cax)
    axis on
    pause(0.01)
end

figure(4)
hold on
se=s{1};
plot(t,double(se(1,:,40)),'r')
