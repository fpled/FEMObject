
r = 10;
L = 1;
l = 1;
H = 2;
r1 = r;
r2 = r*ceil(l/L);
r3 = r*ceil(H/L);

MAT = ELAS_ISOT('E',1,'NU',0.3,'RHO',1);

P1 = POINT([0,0,0]);
P2 = POINT([L,0,0]);
P3 = POINT([L,l,0]);
P4 = POINT([0,0,H]);
P5 = POINT([L,0,H]);
P6 = POINT([L,l,H]);

casmesh = 2;

switch casmesh
    case 1
        S = mesh(DOMAIN(3,P1,P6),r1,r2,r3,'material',MAT);
    case 2
        S = gmsh(DOMAIN(3,P1,P6),L/r);
end
S = setmaterial(S,MAT);
S = final(S);

PL1 = PLAN(P1,P2,P3);
PL2 = PLAN(P4,P5,P6);

S = addcl(S,PL1,'U',0);

K = calc_rigi(S);

%% modes de vibrations
m = 6;
M = calc_mass(S);
[V,D] = calc_mode(K,M,1:m);

%% affichage des modes
ampl = 1;
n = 20;
ampl = cos(linspace(0,2*pi,n))/10*ampl;
ms = [2];
axis0 = [-L/3,L+L/3,-H/3,H+H/3];
nc = ceil(sqrt(length(ms)));
nl = ceil(length(ms)/nc);
figure(1)
for i=1:length(ampl) ;
    for m=1:length(ms)
        if length(ms)>1    
            subplot(nl,nc,m,'replace')  
        else
            clf
        end
        plot(S+ampl(i)*V(:,m),'facecolor','w')
        camlight left
        lighting gouraud
        axis(axis0)
        pause(0.001);
    end
end

pause

%% resolution statique
f = surfload(S,PL2,'FY',1);
% f = zeros(S.nbddl,1);
% [e,numnode] = intersect(S,PL2);
% f(3*numnode) = 0.1*(L*l)/(r1*r2);
% f = SPACEFUN(f,'dual');

q = K\f;

figure(2)
clf
ampl=0.03; 
plot(S+ampl*q,'facecolor','w')

s = calc_sigma(S,q,'node');

figure(3)
plot(s,S,'compo','SMZZ')
colorbar

figure(4)
clf
plot(S,'facecolor','w')
camlight left
lighting gouraud
hold on
vectorplot(S,'F',f,2)
