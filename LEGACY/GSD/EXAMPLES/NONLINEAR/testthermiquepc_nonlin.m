%% Definition probleme
r=15;
P = POINT([ 0,0 ; 1,0 ; 0,2 ; 1,2 ; 1,1 ; 2,1 ; 2,2 ; 0,1 ]);
S1 = mesh(DOMAIN(2,P(1),P(5)),r,r);
S2 = mesh(DOMAIN(2,P(8),P(4)),r,r);
S3 = mesh(DOMAIN(2,P(5),P(7)),r,r);

S = union(S1,S2);
S = union(S,S3);
S = concatgroupelem(S);
RV = RANDVARS();
RV{1}=RVUNIFORM(0.7,1.3);
RV{2}=RVUNIFORM(0.7,1.3);
RV{3}=RVNORMAL(1,0.2);

X = PCMODEL(RV,'order',3,'pcg');

mat=FOUR_ISOT('k',X{1},'k2',X{2});
S = setmaterial(S,mat);
S=final(S);
S=addcl(S,create_boundary(S),'T',0);
f=bodyload(S,[],'QN',10).*X{3};

K0=calc_rigi(S,X);
qlin=pcg(K0,f);
plot(FENODEFIELD(random(qlin)),S)
%% Definition du solveur Newton
N = NEWTONSOLVER('type','modified','tol',1e-7,'tolreact',1e-1,'increment',true);



%% affichage des quantiles
plotcentered = 0;
if plotcentered==1
    qplotquan = qplotquan-expect(qplotquan);
    plotmean=0;
else
    plotmean=1;
end
optionsplot={'surface'};
Ymean = mean(qplotquan);
Y = quantile(qplotquan,[0.1,0.9],1000);
cax = [min(min(Y)),max(max(Y))];

figure(5)
clf
subplot(1,2+plotmean,1)
title('quantile 0.1')
plot(Y{1},S,optionsplot{:})
caxis(cax)
axis off
if ischarin('surface',optionsplot);
    colorbar('horiz')
    zlim(cax)
else
    colorbar
end

subplot(1,2+plotmean,2+plotmean)
title('quantile 0.9')
plot(Y{2},S,optionsplot{:})
caxis(cax)
axis off
if ischarin('surface',optionsplot);
    colorbar('horiz')
    zlim(cax)
else
    colorbar
end
if plotmean
    subplot(1,3,2)
    title('Solution moyenne')
    plot(Ymean,S,optionsplot{:})
    caxis(cax)
    axis off
    if ischarin('surface',optionsplot);
        colorbar('horiz')
        zlim(cax)
    else
        colorbar
    end
end


%% Calcul solution reference non-lineaire
clock0=clock;
qpc = solve(N,f,@(u) calc_fintpc(S,u,X),@(u) calc_rigitangpc(S,u,X),[],[],@(A,b) cgs(A,b,getparam(N,'tol')/100));
etime(clock,clock0)
qradref = spectral_decomposition(qpc,'display');

%% Affichage des modes de la solution non-lineaire : solution reference
figure(13)
clf
plotVwithname(qradref,1:8,'U',S,'surface');


%% Calcul de la solution lineaire par GSD
G = GSDSOLVER('tol',1e-8,'nbfoncmax',10,...
    'display',true,'errorindicator','residual');

[qradlin,result] = solve(G,K0,f);
figure(14)
clf
title('Mode de la solution lineaire')
plotV(qradlin,[],S);


%% Calcul de la solution non-lineaire par GSD : avec reutilisation

G = GSDSOLVER('tol',1e-5,'nbfoncmax',30,'tolini',1e-5,...
    'display',true,'finalSD',true,'righthandSD',true,...
    'errorindicator','residual');
N = setparam(N,'increment',false)
%[qrad,result] = solve(N,f,@(u) calc_fintpc(S,expand(u),X), ...
%    @(u) calc_rigitangpc(S,expand(u),X),[],[],...
%   @(A,b,u) solve(G,A,b,u));
%% SANS REUSE , SANS UPDATE
G = setparam(G,'update',false,'reuse',false);
[qradnoreuse,resultnoreuse] = newton(G,N,f,@(u) calc_fintpc(S,expand(u),X),...
    @(u) calc_rigitangpc(S,expand(u),X),qradlin);
%% SANS REUSE , AVEC UPDATE
G = setparam(G,'update',true,'reuse',false);
[qradupnoreuse,resultupnoreuse] = newton(G,N,f,@(u) calc_fintpc(S,expand(u),X),...
    @(u) calc_rigitangpc(S,expand(u),X),qradlin);

figure(15)
plotaddedfunctions({resultnoreuse,resultupnoreuse})

%% AVEC REUSE , SANS UPDATE
G = setparam(G,'reuse',true,'update',false);
[qradreuse,resultreuse] = newton(G,N,f,@(u) calc_fintpc(S,expand(u),X),...
    @(u) calc_rigitangpc(S,expand(u),X),qradlin);
%% AVEC REUSE , AVEC UPDATE
G = setparam(G,'reuse',true,'update',true);
[qradupreuse,resultupreuse] = newton(G,N,f,@(u) calc_fintpc(S,expand(u),X),...
    @(u) calc_rigitangpc(S,expand(u),X),qradlin);
figure(17)
plotaddedfunctions({resultupreuse,resultupreuse})



%%
GA=setparam(G,'type','arnoldi');
GA = setparam(GA,'reuse',true,'finalSD',true);
GA = setparam(GA,'nbfoncmaxsimul',5,'restart',3,'nbfoncmax',13)

[qradareuse,resultareuse] = newton(GA,N,f,@(u) calc_fintpc(S,expand(u),X),...
    @(u) calc_rigitangpc(S,expand(u),X),qlin);
figure(17)
plotaddedfunctions(resultareuse)

%% Affichage des fonctions ajoutees pendant la resolution GSD incremetale
figure(18)
clf
[addedfun,hh]= plotaddedfunctions({resultnoreuse,resultreuse},1.5);
set(hh(1),'FaceColor','b')
set(hh(2),'FaceColor','y')
legend('No Reuse','Reuse')
xlabel('Newton iteration')
saveas(gcf,'E:\REDACTION\matlab_results\thermique_nonlin_addedfun_reuse_noreuse.eps','epsc2')


%% Affichage des modes de la solution non-lineaire : solution GSD

figure(15)
clf
title('Mode de la solution non lineaire : newton et solver lineaire GSD')
optionsplot = {'surface'};
plotV(qrad,[],S,optionsplot{:});
if ischarin('surface',optionsplot)
    saveas(gcf,'E:\REDACTION\matlab_results\thermique_nonlin_P-GSD_newton_spacefunctions_surf.eps','epsc2')
else
    saveas(gcf,'E:\REDACTION\matlab_results\thermique_nonlin_P-GSD_newton_spacefunctions.eps','epsc2')
end

%% Calcul de la solution non-lineaire par GSD : Newton global

NG = NEWTONSOLVER('type','modified','tol',1e-2,...
    'tolreact',1e-1,'increment',true,'display',false);
G = GSDSOLVER('tol',1e-4,'nbfoncmax',8,'pfixtol',5e-2,...
    'display',true,'update',false,'orthocrit',1e-12);
G = setparam(G,'type','power');
[qradng,resultng] = newtonglob(G,NG,S,f,X,'reference',qpc);

%%
figure(10021)
semilogy(resultng.error,'k^-','LineWidth',2,'MarkerSize',10)
fsize=18;
xlim([0.5,length(resultng.error)+0.5])
set(gca,'FontSize',18)
xlabel('order M')
ylabel('error')
saveas(gcf,'E:\REDACTION\matlab_results\thermique_nonlin_PU-GSD_newtonglobal_error.eps','epsc2')

%%
G = setparam(G,'type','arnoldi');
G = setparam(G,'inittype','allone');
G = setparam(G,'restart',0);
G = setparam(G,'nbfoncmax',4);
[qradng,resultng_a] = newtonglob(G,NG,S,f,X,'reference',qpc);

figure(10021)
hold on
semilogy(resultng_a.error,'rs-','LineWidth',2,'MarkerSize',10)
legend('PU-GSD','AGSD')
%% Affichage des modes de la solution non-lineaire par GSD : Newton global
figure(16)
clf
title('Mode de la solution non lineaire : solveur GSD non-lineaire')
U = qradng;
plotV(U,[],S,'surface');
colormap('default')
saveas(gcf,'E:\REDACTION\matlab_results\thermique_nonlin_A-GSD_newtonglobal_spacefunctions_surf.eps','epsc2')

figure(160)
clf
title('Mode de la solution non lineaire : solveur GSD non-lineaire')
plotV(U,[],S);
%colormap(flipud(gray))
saveas(gcf,'E:\REDACTION\matlab_results\thermique_nonlin_A-GSD_newtonglobal_spacefunctions.eps','epsc2')


%% affichage des quantile de la temperature esur une droite
qplotquan = qradng;
%qplotquan = qradreuse-expect(qradreuse);
qplotquan = orthogonalizeL(qplotquan);
qplotquan=normsort(qplotquan);
choixdroite=1;
switch choixdroite
case 1
    D10 = LIGNE(POINT([0.5,0]),POINT([0.5,2]));
    [rep,P10]=ispointin(D10,POINT(S.node));
    t=double(getcoord(P10));t=t(1,2,:);t=t(:);
case 2
    L10 = LIGNE(POINT([0.5,0]),POINT([0.5,1.5]));
    L11 = LIGNE(POINT([0.5+1/r,1.5]),POINT([2,1.5]));
    [rep,P10]=ispointin(L10,POINT(S.node));
    [rep,P11]=ispointin(L11,POINT(S.node));
    P10 = [P10;P11];
    t=[0:1/r:1.5*2];t=t(:);

end
qdroitequan=eval_sol(S,qplotquan,P10,'T');


m = getm(qplotquan);
nl = ceil(sqrt(m));nc = ceil(m/nl);

Ymean = expect(qdroitequan);
fillcolor={'k','b','r','m','g','c','y'};
fillcolor=fillcolor(end:-1:1);
nr = 10000;
quan = [0.05,0.95];
figure(200*choixdroite)
clf
scanmode=[4:-1:1];
nl = ceil(sqrt(length(scanmode)));nc = ceil(length(scanmode)/nl);
for k=1:length(scanmode)
    i=scanmode(k);
    V = getmode(qdroitequan,i);
    Y = quantile(V,quan,nr);
    subplot(nl,nc,i)
    plotenveloppe(t,Y);
    title(['mode ' num2str(i)])
end

figure(201*choixdroite)
clf
hold on
leg = {};
scanmode = [4:-1:1];
for k=1:length(scanmode)
    kc = mod(k-1,length(fillcolor))+1;
    i = scanmode(k);
    V = getmode(qdroitequan,i);
    Y = quantile(V,quan,nr);
    plotenveloppe(t,Y,fillcolor{kc});
    leg = [leg , {['mode ' num2str(i)]}];
end
leg = [leg,'mean solution'];
plot(t,Ymean,'k-','LineWidth',2);
legend(leg{:})
title('envelep of modes')

figure(202*choixdroite)
clf
hold on
leg = {};
scanmode = [6,3:-1:1];
for k=1:length(scanmode)
    kc = mod(k-1,length(fillcolor))+1;
    i=scanmode(k);
    V = getmode(qdroitequan,1:i);
    Y = quantile(V,quan,nr);
    plotenveloppe(t,Y,fillcolor{kc});
    leg = [leg , {['mode 1 to ' num2str(i)]}];
end
plot(t,Ymean,'k-','LineWidth',2)
leg = [leg,'mean solution'];
legend(leg{:})

