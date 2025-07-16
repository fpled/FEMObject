S = gmsh2femobject(3,'gmsh_biplates.msh',2);
S = convertelem(S,'DKT');
mat = ELAS_SHELL('E',1,'NU',.3,'RHO',1,'DIM3',.1);
S = setmaterial(S,mat);
S = final(S);

figure(1)
clf
plot(S)

figure(2)
clf
plotfacets(S);

figure(3)
clf
plotridges(S)

nodeefforts = getnumber(getnode(getfacet(S,6)));
nodedirichlet = getnumber(getnode(getfacet(S,4)));

nodeoutput = getnumber(getnode(getridge(S,3)));
ddloutput = findddl(S,'UY',nodeoutput,'free');

figure(4)
clf
plot(S)
plot(getnode(S,nodeefforts),'.r')
plot(getnode(S,nodedirichlet),'.g')
plot(getnode(S,nodeoutput),'.b')

figure(5)
clf
plotparamelem(S,'group');

S = addcl(S,nodedirichlet,'U');
f = nodalload(S,nodeefforts,'FY',1);

K1 = calc_rigi(S,'selgroup',1); % matrice de rigidite du groupe d'element 1
K2 = calc_rigi(S,'selgroup',2); % matrice de rigidite du groupe d'element 2
K = calc_rigi(S); % matrice de rigidite globale K = K1+K2
M = calc_mass(S);

q=K\f;

figure(6)
clf
plot_sol(S,q,'displ',1,'solid')

figure(7)
clf
plot_sol(S,q,'displ',2,'solid')

figure(8)
clf
plot_sol(S,q,'displ',3,'solid')

figure(9)
clf
plot_sol(S,q,'rotation',1,'solid')

figure(10)
clf
plot_sol(S,q,'rotation',2,'solid')

figure(11)
clf
plot_sol(S,q,'rotation',3,'solid')

%% calcul vibrations sans amort
[V,D] = calc_mode(K,M,1:4);

figure(12)
clf
for i=1:size(V,2)
    ampl = 1/max(abs(V(:,i)));
    subplot(2,2,i)
    plot(S,'Color','none','FaceColor','k','FaceAlpha',0.5)
    plot(S+ampl*V(:,i),'Color','none','FaceColor','r','FaceAlpha',0.5)
    lighting gouraud
end

%%
for i=1:3
    E1 = 1;
    E2 = random(RVUNIFORM(1,2));
    K = E1*K1+E2*K2;
    C = .01*K+.01*M;
    scanw = linspace(0,2,10);
    [urw,uiw,umod,uphi] = calc_frf(scanw,K,C,M,f,ddloutput);
    
    figure(13)
    hold on
    plot(scanw,log(umod),'Color',getfacecolor(i+1));

    w = 1.7;
    u = (K-w^2*M+j*C*w)\f;
    
    figure(14)
    hold on
    plot(S+.001*real(u),'FaceColor',getfacecolor(i+1),'FaceAlpha',0.5)
    
    figure(15)
    hold on
    plot(S+.001*imag(u),'FaceColor',getfacecolor(i+1),'FaceAlpha',0.5)
end
