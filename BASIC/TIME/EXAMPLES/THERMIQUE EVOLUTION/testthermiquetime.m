

r1=15;r2=15;
nt=40;k=46;
T = 1e5;
t = linspace(0,T,nt+1);
dt=T/nt ;

close all
display_=1;
L=1;H=1;
%k=46;
rho=7800;
C=460;
mat=FOUR_ISOT('c',rho*C,'k',k);

P1=POINT([0,0]);
P2=POINT([L,0]);
P3=POINT([L,H]);
P4=POINT([0,H]);
P5=POINT([0,H/2]);
P6 = POINT([L/2,0]);
S=mesh(DOMAIN(2,P1,P3),r1,r2,mat);
%S=setmaterial(S,mat);
S=final(S);


D1  = STRAIGHTLINE(P1,P4);
D2  = STRAIGHTLINE(P2,P3);
D3  = STRAIGHTLINE(P1,P2);
D4  = STRAIGHTLINE(P4,P3);
L1  = LINE(P1,P5);
L2  = LINE(P1,P6);

%S=addcl(S,D2,'T',0);
S=addcl(S,D4,'T',0);
%S=addcl(S,D3,'T',0);


K=calc_rigi(S,mat);
M=calc_mass(S,mat);

% fun1 = inline('x(:,2)','x');
% fun2 = inline('-x(:,2)','x');
f=surfload(S,L1,{'QN'},1e4);
f=surfload(S,L2,{'QN'},1e4);

q=K\f;

q=unfreevector(S,q);


if display_

    figure(1)
    clf

    plot(FENODEFIELD(q),S)
    title('solution stationnaire')
    axis image
    cax0=caxis;
end


R = chol(M/dt+K);

f = f*ones(1,nt);
disp('resolution temporelle')
qt{1}=zeros(size(K,1),1);
for i=1:nt
    pourcentage(i,nt)
    qt{i} = R\(R'\(f(:,i)+1/dt*M*qt{i-(i>1)}*(i>1)) );   
end


if display_
    figure(2)
    for i=1:5:nt
        clf   
        plot(FENODEFIELD(unfreevector(S,qt{i})),S);
        title(['time ' num2str(t(i)) ' s'])
        caxis(cax0)
        axis image
        colorbar
        pause(0.1)

    end
end
