
% S = gmsh2femobject(2,[getfemobjectoptions('path') 'APPLIS/TUTORIALS/evolution_problems/example_canister.geo'],2);
cl1 = 0.02;
cl2 = 0.04;
cl3 = 0.02;
cltip = 0.01;
S = gmshCanister(cl1,cl2,cl3,cltip,'gmsh_canister');

figure(1)
clf
plot(S)

figure(2)
clf
plotfacets(S)

Sc = S;


%% probleme et maillage
figure(44)
clf;plot(create_boundary(S));
h1=plot(S,'selgroup',1,'FaceColor','c','EdgeColor','none');
h2=plot(S,'selgroup',2:3,'FaceColor','y','EdgeColor','none');
set(gca,'FontSize',16)
legend([h1(1),h2(2)],'\Omega_1','\Omega_2')
% myprint('canister','domain','jpeg');


figure(45)
clf;plot(S);
h1=plot(S,'selgroup',1,'FaceColor','c','EdgeColor','none');
h2=plot(S,'selgroup',2:3,'FaceColor','y','EdgeColor','none');
set(gca,'FontSize',12)
legend([h1(1),h2(1)],'\Omega_1','\Omega_2')
% myprint('canister','mesh','jpeg');

%%
figure(46)
clf;plot(create_boundary(keepgroupelem(S,1)));
plot(create_boundary(keepgroupelem(S,2:3)));
YY=18;
text(0,0,'\Gamma_1','FontSize',YY)
text(0,0,'\Gamma_2','FontSize',YY)
text(0,0,'\Omega_1','FontSize',YY)
text(0,0,'\Omega_2','FontSize',YY)

% myprint('canister','problem','jpeg');
ylim([0,1.8])

%%

mat = FOUR_ISOT('k',1);
Sc = setmaterial(Sc,mat);
Sc = final(Sc);
Sc = addcl(Sc,POINT([0,0]),'T');
P5 = POINT(getnode(getridge(Sc,5)));
P6 = POINT(getnode(getridge(Sc,6)));
P15 = POINT(getnode(getridge(Sc,15)));
P16 = POINT(getnode(getridge(Sc,16)));
L1 = LIGNE(P5,P6);
L2 = LIGNE(P15,P16);
f1 = surfload(Sc,L1,'QN',-1);
f2 = surfload(Sc,L2,'QN',1);
f = f1+f2;
K = calc_rigi(Sc);
phi = K\f;
v = FENODEFIELD(calc_sigma(Sc,phi,'node'));
V = getvalue(v);
V = {{FENODEFIELD(V(:,1)),FENODEFIELD(V(:,2))}};
figure(4)
clf
plot(phi,Sc)
hold on
quiver(v,Sc,6,'k')
ylim([0,1.7])
% myprint('C:/','example_canister_flow',{'jpeg'});

%%
mat = FOUR_ISOT('k',1,'b',V,'c',1);
S = setmaterial(S,mat);
S = final(S);
S = addcl(S,[],'T');

[L1mesh,node1] = intersect(S,L1);
u0 = zeros(getnbddl(S),1);
u0(node1) = 1;

kreac1 = 0.1;
kreac2 = 10;
kdiff = 0.01;
kadv = 2;
ku0 = 1;

Mx = calc_matrix(S,@mass);
Adiff = calc_matrix(S,@diff);
Aadv = calc_matrix(S,@adv,'intorder','mass');
Areac1 = calc_matrix(S,@mass,'selgroup',1);
Areac2 = calc_matrix(S,@mass,'selgroup',2:3);

foutput = bodyload(keepgroupelem(S,2),[],'QN',1,'free');

timespeed = 0;

t0 = 0;
nt = 100;
if timespeed == 0
    t1 = 2;
    T = TIMEMODEL(t0,t1,nt);
    N = DGTIMESOLVER(T,1);
    A = kdiff*Adiff+kadv*Aadv+kreac1*Areac1+kreac2*Areac2;
    b0 = -A*(ku0*u0);
    A = freematrix(S,A);
    b0 = freevector(S,b0);
    b = b0*one(N);

    u = A\b0;
    u = unfreevector(S,u)+ku0*u0;

    figure(6)
    clf
    plot(u,S,'surface');

else
    t1 = 4;
    T = TIMEMODEL(t0,t1,nt);
    N = DGTIMESOLVER(T,1);
    A = kdiff*Adiff+kreac1*Areac1+kreac2*Areac2;
    A = A + N(@(t) 2- 4*t/gett1(N).*(1-t/gett1(N)))*(kadv*Aadv); %(1+cos(N,pi/gett1(N)))
    b = -A*(ku0*u0);
    A = freematrix(S,A);
    b = freevector(S,b);
end

Mx = freematrix(S,Mx);
%

N = setparam(N,'display',true);
uref = dsolve(N,b,Mx,A);

utplot = unfreevector(S,uref)+ku0*u0;

figure(8)
clf
evol(utplot,S,'rescale','z')

figure(9)
clf
plot(foutput'*uref);
% ftoutput = one(N);
ftoutput = getfinalequivalent(N);
fxoutput = foutput;

if timespeed
    fichname = 'example_canister_timespeed_';
else
    fichname = 'example_canister_';
end
save(fichname,'A','b','Mx','uref','N','S')
