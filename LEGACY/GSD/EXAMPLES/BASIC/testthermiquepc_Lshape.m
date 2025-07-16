r=10;
P = POINT([ 0,0 ; 1,0 ; 0,2 ; 1,2 ; 1,1 ; 2,1 ; 2,2 ]);
S1 = mesh(DOMAIN(2,P(1),P(4)),r,r*2);
S2 = mesh(DOMAIN(2,P(5),P(7)),r,r);
S = union(S1,S2);
S = concatgroupelem(S);

RV=RANDVARS();

secdet = 1;

e=0.5;
%kmarg = RFMARGINAL(RVGAMMA(1/0.25^2,0.25^2));
kmarg = RFMARGINAL(RVUNIFORM(0.5,1.5));
k = RANDFIELD(kmarg,EXP2CORREL(2));
m=3;p=4;msd=10;
if m>0
    kmat = KL(k,m,S,'order',p);
    kmat = spectral_decomposition(kmat,'nbfoncmax',msd);
else
    kmat = mean(getrandvar(kmarg));
end

if secdet
    RV = RANDVARS(RANDVARS(kmat));
else
    RV = RANDVARS(RANDVARS(kmat),RVNORMAL(0.5,0.2),RVNORMAL(0,0.2));
end
X = PCMODEL(RV,'order',p,RANDPOLYS(RV));
if m>0
    kmat = project(kmat,X);
end



L1 = LIGNE(P(6),P(7));
L2 = LIGNE(P(1),P(3));
L3 = LIGNE(P(3),P(7));
L4 = LIGNE(P(6),P(5));
L5 = LIGNE(P(2),P(5));
L6 = LIGNE(P(1),P(2));

mat=setnumber(FOUR_ISOT('k',kmat),1);
S = setmaterial(S,mat);
S=final(S);
S=addcl(S,L1,'T',0);
S=addcl(S,L2,'T',0);
S=addcl(S,L3,'T',0);
S=addcl(S,L4,'T',0);
S=addcl(S,L5,'T',0);

Ksto = calc_rigipc(S,X);


f1 = bodyload(S,[],'QN',1);
f2 = surfload(S,L6,'QN',1);
if secdet
    fsto = f1+f2;
else
    fsto = f1*X{end-1}+f2.*X{end};
end

disp('Resolution pcg')
tic
if israndom(Ksto)
    qpc = pcg(Ksto,fsto,1e-6);
else
    qpc = Ksto\fsto;
end
toc
time_pcg=toc;



%% Resolution spectrale
GSD = GSDSOLVER('tol',1e-4,'nbfoncmax',6,...
    'display',true,'update',true,'errorindicator','residual','type','arnoldi','restart',3);
[qrad,result]=solve(GSD,Ksto,fsto);
%%
GSD = GSDSOLVER('tol',1e-4,'nbfoncmax',6,...
    'display',true,'update',true,'errorindicator','residual','type','power');
[qrad,result]=solve(GSD,Ksto,fsto);

%%
GSD = GSDSOLVER('tol',1e-4,'pfixtol',1e-8,'pfixmax',4,'nbfoncmax',4,...
    'display',true,'update',true,...
    'errorindicator','residual','type','power',...
    'subspaceiteration',0,'direct',true);
[qrad,result]=solve(GSD,Ksto,fsto);

%% subspace methods
GSD = GSDSOLVER('tol',1e-4,'pfixtol',1e-8,'pfixmax',5,'nbfoncmax',3,...
    'display',true,'update',true,...
    'errorindicator','residual','type','powersubspace',...
    'direct',true);
[qrad,result]=solve(GSD,Ksto,fsto);


%%


clear result_p
clear result_ps
clear leg

fsizetransp = 24;
fsize=16;

courbestyle = {'b*-','rs-','ko-','m+-','gd-','b^-','r<-','k>-','mx-'};
leg={};
figure(11)
clf
M=6;

GSD = GSDSOLVER('tol',1e-4,'pfixtol',1e-3,'pfixmax',3,'nbfoncmax',M,...
    'display',true,'update',true,...
    'errorindicator','rayleigh','type','power',...
    'subspaceiteration',0,'direct',false);

GSD=setparam(GSD,'update',false);

scansub=[0:2];
for k=1:length(scansub)
    GSD=setparam(GSD,'subspaceiteration',scansub(k));
    for i=1:M
        GSD=setparam(GSD,'nbfoncmax',i);
        [qrad,result_ps{k}{i}]=solve(GSD,Ksto,fsto,[],'reference',qpc);
        error_ps{k}(i) = result_ps{k}{i}.error(end);
    end
    semilogy(error_ps{k},courbestyle{length(leg)+1},'LineWidth',2,'MarkerSize',10)
    hold on
    leg = [leg , {['PS(power)-GSD, iter='  num2str(scansub(k)) ]}];
    pause(0.01)
end

%%
GSD=setparam(GSD,'update',true);
clear result_psu
clear error_psu
scansub=[0:2];
for k=1:length(scansub)
    GSD=setparam(GSD,'subspaceiteration',scansub(k));
    for i=1:M
        GSD=setparam(GSD,'nbfoncmax',i);
        [qrad,result_psu{k}{i}]=solve(GSD,Ksto,fsto,[],'reference',qpc);
        error_psu{k}(i) = result_psu{k}{i}.error(end);
    end
    semilogy(error_psu{k},courbestyle{length(leg)+1},'LineWidth',2,'MarkerSize',10)
    hold on
    leg = [leg , {['PS(power-update)-GSD, iter='  num2str(scansub(k)) ]}];
    pause(0.01)
end


h=gca;
set(gca,'Fontsize',fsize)
xlabel('order M','FontSize',fsize)
ylabel('error','FontSize',fsize)
legend(leg{:})

%%

GSD = GSDSOLVER('tol',1e-4,'nbfoncmax',i,...
    'display',true,...
    'errorindicator','residual','type','arnoldi',...
    'subspaceiteration',1,'direct',true);
clear result_a
clear error_a
scansub=[0:2];
M=6;
for k=1:length(scansub)
    for i=1:M
        GSD=setparam(GSD,'nbfoncmax',i);
        GSD=setparam(GSD,'subspaceiteration',scansub(k));
        [qrad,result_a{k}{i}]=solve(GSD,Ksto,fsto,[],'reference',qpc);
        error_a{k}(i)=result_a{k}{i}.error(end);
    end
    semilogy(error_a{k},courbestyle{length(leg)+1},'LineWidth',2,'MarkerSize',10)
    leg = [leg , {['PS(arnoldi)-GSD, iter=' num2str(scansub(k))]}];
end
legend(leg{:})
%% POWER SUBSPACE
figure(31)
clf
leg={};

GSD = GSDSOLVER('tol',1e-4,'pfixtol',1e-3,'pfixmax',1,'nbfoncmax',6,...
    'display',true,...
    'errorindicator','rayleigh','type','powersubspace',...
    'direct',true);

clear result_pss
clear error_pss
scansub=[0:2];
M=6;
for k=1:length(scansub)
    for i=1:M
        GSD=setparam(GSD,'nbfoncmax',i);
        GSD=setparam(GSD,'pfixmax',scansub(k));
        [qrad,result_pss{k}{i}]=solve(GSD,Ksto,fsto,[],'reference',qpc);
        error_pss{k}(i)=result_pss{k}{i}.error(end);
    end
    semilogy(error_pss{k},courbestyle{length(leg)+1},'LineWidth',2,'MarkerSize',10)
    hold on
    leg = [leg , {['PS-GSD, iter=' num2str(scansub(k))]}];
end
legend(leg{:})


%%



test=5
switch test
    case 1
        A = Ksto ;
        b = fsto;
        qref = qpc ;
    case 2
        A = expect(Ksto);
        b = A*qrad;
        qref = qrad;
    case 3
        b = expect(fsto);
        A = Ksto;
        qref = pcg(A,b,1e-10);
    case 4
        b = qrad;
        A = Ksto;
        qref = pcg(A,b,1e-10);
    case 5
        A=expect(Ksto);
        Asup = triu(A,1);Ainf=tril(A,-1);Adiag = diag(diag(A));
        A = Asup/2*0+Ainf+Adiag;
        b = A*qrad ;
        qref=qrad;
end


%%
GSD = GSDSOLVER('tol',1e-4,'nbfoncmax',6,...
    'display',false,'pfixmax',10,'pfixtol',1e-7,...
    'errorindicator','reference','type','power','direct',true,...
    'update',true,'subspaceiteration',5);
[qradsol,result]=solve(GSD,A,b,[],'reference',qref);
result.error

%%
GSD = setparam(GSD,'subspaceiteration',0);
GSD = setparam(GSD,'update',true);
[qrad1,result1]=solve(GSD,A,b,[],'reference',qref);
result1.error
GSD = setparam(GSD,'update',false);
[qrad2,result2]=solve(GSD,A,b,[],'reference',qref);
result2.error


%% Resolution spectrale
for i=0:3
    GSD = GSDSOLVER('tol',1e-4,'pfixtol',1e-10,'pfixmax',15,'nbfoncmax',4,...
        'display',true,'update',false,...
        'errorindicator','residual','type','power',...
        'subspaceiteration',i,'direct',true);
    [qrad,result]=solve(GSD,Ksto,fsto,[],'reference',qpc);
    err(i+1) = norm(qpc-qrad,Ksto)/norm(qpc,Ksto)
    errerr(i+1) = (sqrt(norm(qpc,Ksto)^2-result.rayg{end})-norm(qpc-qrad,Ksto))/norm(qpc-qrad,Ksto)
end
figure(434)
semilogy(err,'*-')
title('error')
figure(435)
semilogy(errerr,'r-^')
title('error on error estimate')


%% Description du probleme : figures
figure(200)
clf
plot(S,'FaceColor',[0.6,0.9,0.9])
%plot(S)
colu = 'b' ;
colg = 'r' ;
linwid = 3;
plot(L1,'LineWidth',linwid,'EdgeColor',colu)
plot(L2,'LineWidth',linwid,'EdgeColor',colu)
plot(L3,'LineWidth',linwid,'EdgeColor',colu)
plot(L4,'LineWidth',linwid,'EdgeColor',colu)
plot(L5,'LineWidth',linwid,'EdgeColor',colu)
plot(L6,'LineWidth',linwid,'EdgeColor',colg)

annotation('textarrow',[0.7 0.65],[0.35 0.5],'String','\partial_1\Omega','FontSize',18)
annotation('textarrow',[0.55 0.45],[0.2 0.12],'String','\partial_2\Omega','FontSize',18)

%saveas(gcf,'E:\REDACTION\matlab_results\Lshape_mesh_boundaryconditions.eps','epsc2')
%saveas(gcf,'E:\REDACTION\matlab_results\Lshape_mesh_boundaryconditions_sanscol.eps','epsc2')

figure(201)
clf
%plot(S,'FaceColor',[0.6,0.9,0.9])
plot(S)
repnode=330;
plot(S.node(repnode),'Marker','bx','MarkerSize',12)
plottext(S.node(repnode),'P_1','Color','b','FontSize',18)
repnode=640;
plot(S.node(repnode),'Marker','bx','MarkerSize',12)
plottext(S.node(repnode),'P_2','FontSize',18,'Color','b')
repnode=1097;
plot(S.node(repnode),'Marker','bx','MarkerSize',12)
plottext(S.node(repnode),'P_3','FontSize',18,'Color','b')


%saveas(gcf,'E:\REDACTION\matlab_results\Lshape_mesh_boundaryconditions.eps','epsc2')
%saveas(gcf,'E:\REDACTION\matlab_results\Lshape_mesh_withpoints.eps','epsc2')

figure(202)
clf
plot(S,'FaceColor',[0.6,0.9,0.9])
%plot(S)
colu = 'b' ;
colg = 'r' ;
linwid = 3;
plot(L1,'LineWidth',linwid,'EdgeColor',colu)
plot(L2,'LineWidth',linwid,'EdgeColor',colu)
plot(L3,'LineWidth',linwid,'EdgeColor',colu)
plot(L4,'LineWidth',linwid,'EdgeColor',colu)
plot(L5,'LineWidth',linwid,'EdgeColor',colu)
plot(L6,'LineWidth',linwid,'EdgeColor',colg)

annotation('textarrow',[0.7 0.65],[0.35 0.5],'String','\partial_1\Omega','FontSize',18)
annotation('textarrow',[0.55 0.45],[0.2 0.12],'String','\partial_2\Omega','FontSize',18)

repnode=330;
plot(S.node(repnode),'Marker','kx','MarkerSize',12)
plottext(S.node(repnode),'P_1','Color','k','FontSize',18)
repnode=640;
plot(S.node(repnode),'Marker','kx','MarkerSize',12)
plottext(S.node(repnode),'P_2','FontSize',18,'Color','k')
repnode=1097;
plot(S.node(repnode),'Marker','kx','MarkerSize',12)
plottext(S.node(repnode),'P_3','FontSize',18,'Color','k')

%saveas(gcf,'E:\REDACTION\matlab_results\Lshape_mesh_boundaryconditions_withpoints.eps','epsc2')
%saveas(gcf,'E:\REDACTION\matlab_results\Lshape_mesh_boundaryconditions_sanscol.eps','epsc2')

%% -------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% CALCUL DES ERREURS FONCTION DU TEMPS
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------

%% erreur pcg en fonction du temps de resolution
errorpcg = [5e-1,10.^([-1:-1:-4])];
for k=1:length(errorpcg)
    tic
    fprintf('erreur = %d',errorpcg(k))
    qpctemp = pcg(Ksto,fsto,errorpcg(k));
    timepcg(k)=toc
    %errorpcg(k) = norm(qpctemp-qpc)/norm(qpc)
end

timepcgp{p}=timepcg;
errorpcgp{p}=errorpcg;
timepcgr{r}=timepcg;
errorpcgr{r}=errorpcg;


%%


GSD = GSDSOLVER('tol',1e-4,'pfixtol',1e-3,'pfixmax',2,'nbfoncmax',10,...
    'display',true,'update',true,...
    'errorindicator','rayleigh','type','power',...
    'subspaceiteration',0,'direct',false);
[qrad,result_powertime]=solve(GSD,Ksto,fsto);
timepower=result_powertime.time;
[qrad,result_power]=solve(GSD,Ksto,fsto,[],'reference',qpc);
errorpower = result_power.error ;

semilogy(timepower,errorpower,'b*-','LineWidth',2)

timepowerp{p}=timepower;
errorpowerp{p}=errorpower;

%%

GSD = GSDSOLVER('tol',1e-4,'nbfoncmax',4,...
    'errorindicator','none','type','arnoldi','direct',false,'display','true',...
    'inittype','one','orthocrit',1e-10);
for i=1:8
    GSD = setparam(GSD,'nbfoncmax',i+1);
    GSD = setparam(GSD,'errorindicator','none');
    [qrada,result_atime]=solve(GSD,Ksto,fsto);
    timea(i)=result_atime.time(end);
    qradasd = spectral_decomposition(qrada,'nbfoncmax',i,'pfixmax',7,'pfixtol',1e-10);
    errora(i)=norm(qradasd-qpc)/norm(qpc);
    
    
end

semilogy(timea,errora,'b*-','LineWidth',2)

timeap{p}=timea;
errorap{p}=errora;
timear{r} = timea;
errorar{r}=errora;


%%

fsizetransp = 20;
fsize=16;

%% comparaison arnoldi power pcg
rr=3;
figure(10)
clf
semilogy(timepcgp{3}(1:rr),errorpcgp{3}(1:rr),'k*--','LineWidth',3,'MarkerSize',10)
hold on
semilogy(timepowerp{3},errorpowerp{3},'rs--','LineWidth',3,'MarkerSize',10)
semilogy(timeap{3},errorap{3},'b>--','LineWidth',3,'MarkerSize',10)
semilogy(timepcgp{5}(1:rr),errorpcgp{5}(1:rr),'k*-','LineWidth',3,'MarkerSize',10)
semilogy(timepowerp{5},errorpowerp{5},'rs-','LineWidth',3,'MarkerSize',10)
semilogy(timeap{5},errorap{5},'b>-','LineWidth',3,'MarkerSize',10)
semilogy(timepcgp{}(1:rr),errorpcgp{7}(1:rr),'k*:','LineWidth',3,'MarkerSize',10)
semilogy(timepowerp{7},errorpowerp{7},'rs:','LineWidth',3,'MarkerSize',10)
semilogy(timeap{7},errorap{7},'b>:','LineWidth',3,'MarkerSize',10)
ylim([10^-4,1])


h=gca;
set(gca,'Fontsize',fsizetransp)
xlabel('time (s)','FontSize',fsizetransp)
ylabel('error','FontSize',fsizetransp)
legend('PCG , p=3','P-GSD , p=3','A-GSD , p=3','PCG , p=5','P-GSD , p=5','A-GSD , p=5','PCG , p=7','P-GSD , p=7','A-GSD , p=7')

saveas(gcf,'E:\REDACTION\matlab_results\Lshaperandomfield\Lshape_errorvstime_arnoldi_power_pcg_orderp.eps','epsc2')

%% comparaison arnoldi pcg
rr=3;
figure(11)
clf
semilogy(timepcgp{3}(1:rr),errorpcgp{3}(1:rr),'k*--','LineWidth',3,'MarkerSize',10)
hold on
semilogy(timeap{3},errorap{3},'b>--','LineWidth',3,'MarkerSize',10)
semilogy(timepcgp{5}(1:rr),errorpcgp{5}(1:rr),'k*-','LineWidth',3,'MarkerSize',10)
semilogy(timeap{5},errorap{5},'b>-','LineWidth',3,'MarkerSize',10)
semilogy(timepcgp{7}(1:rr),errorpcgp{7}(1:rr),'k*:','LineWidth',3,'MarkerSize',10)
semilogy(timeap{7},errorap{7},'b>:','LineWidth',3,'MarkerSize',10)
ylim([10^-4,1])
h=gca;
set(gca,'Fontsize',fsizetransp)
xlabel('time (s)','FontSize',fsizetransp)
ylabel('error','FontSize',fsizetransp)
legend('PCG , p=3','A-GSD , p=3','PCG , p=5','A-GSD , p=5','PCG , p=7','A-GSD , p=7')

saveas(gcf,'E:\REDACTION\matlab_results\Lshaperandomfield\Lshape_errorvstime_arnoldi_pcg_orderp.eps','epsc2')


%% comparaison arnoldi pcg
rr=4;
scanr=[10,20,30];
scann=[341,1282,2822];
figure(1111)
clf
courbestyle1={'k*--','k*-','k*:'};
courbestyle2={'b>--','b>-','b>:'};


leg={};
for kkk=1:length(scanr)
    semilogy(timepcgr{scanr(kkk)}(1:rr),errorpcgr{scanr(kkk)}(1:rr),courbestyle1{kkk},'LineWidth',3,'MarkerSize',10)
    leg=[leg ,{['PCG n=' num2str(scann(kkk))]}];
    hold on
    semilogy(timear{scanr(kkk)},errorar{scanr(kkk)},courbestyle2{kkk},'LineWidth',3,'MarkerSize',10)
    leg=[leg ,{['A-GSD n=' num2str(scann(kkk))]}];
end
ylim([10^-4,1])
h=gca;
set(gca,'Fontsize',fsizetransp)
xlabel('time (s)','FontSize',fsizetransp)
ylabel('error','FontSize',fsizetransp)
legend(leg{:})

saveas(gcf,'E:\REDACTION\matlab_results\Lshaperandomfield\Lshape_errorvstime_arnoldi_pcg_sizen.eps','epsc2')


%% comparaison arnoldi pcg
rr=3;
figure(101)
clf
semilogy(timepcgp{5}(1:rr),errorpcgp{5}(1:rr),'k*-','LineWidth',3,'MarkerSize',10)
hold on
semilogy(timeap{5},errorap{5},'b>-','LineWidth',3,'MarkerSize',10)
ylim([10^-4,1])
xlim([0,50])
h=gca;
set(gca,'Fontsize',fsizetransp)
xlabel('time (s)','FontSize',fsizetransp)
ylabel('error','FontSize',fsizetransp)
legend('PCG','A-GSD')

saveas(gcf,'E:\REDACTION\matlab_results\Lshaperandomfield\Lshape_errorvstime_arnoldi_pcg_orderp.eps','epsc2')

%% comparaison power pcg
rr=3;
figure(12)
clf
semilogy(timepcgp{3}(1:rr),errorpcgp{3}(1:rr),'k*--','LineWidth',3,'MarkerSize',10)
hold on
semilogy(timepowerp{3},errorpowerp{3},'rs--','LineWidth',3,'MarkerSize',10)
semilogy(timepcgp{5}(1:rr),errorpcgp{5}(1:rr),'k*-','LineWidth',3,'MarkerSize',10)
semilogy(timepowerp{5},errorpowerp{5},'rs-','LineWidth',3,'MarkerSize',10)
semilogy(timepcgp{7}(1:rr),errorpcgp{7}(1:rr),'k*:','LineWidth',3,'MarkerSize',10)
semilogy(timepowerp{7},errorpowerp{7},'rs:','LineWidth',3,'MarkerSize',10)
ylim([10^-4,1])


h=gca;
set(gca,'Fontsize',fsizetransp)
xlabel('time (s)','FontSize',fsizetransp)
ylabel('error','FontSize',fsizetransp)
legend('PCG , p=3','P-GSD , p=3','PCG , p=5','P-GSD , p=5','PCG , p=7','P-GSD , p=7')
saveas(gcf,'E:\REDACTION\matlab_results\Lshaperandomfield\Lshape_errorvstime_power_pcg._orderp.eps','epsc2')


% pour charger les temps : load('Lshape_errorvstime_n20')


%%

qpcsd = spectral_decomposition(qpc,'nbfoncmax',8,'pfixmax',7,'pfixtol',1e-10);
%%
for i=1:getm(qpcsd)
    errorqpcsd(i)=norm(truncate(qpcsd,1:i)-qpc)/norm(qpc);
    errorAqpcsd(i)=norm(truncate(qpcsd,1:i)-qpc,Ksto)/norm(qpc,Ksto);
end


%% POWER SUBSPACE

GSD = GSDSOLVER('tol',1e-4,'pfixtol',1e-14,'pfixmax',11,'nbfoncmax',8,...
    'display',true,...
    'errorindicator','rayleigh','type','powersubspace',...
    'direct',true);
for i=1:8
    GSD = setparam(GSD,'nbfoncmax',i);
    [qrad,resultps]=solve(GSD,Ksto,fsto);
    errorpsgsd(i)=norm(qrad-qpc)/norm(qpc);
    errorApsgsd(i)=norm(qrad-qpc,Ksto)/norm(qpc,Ksto);
end
%%
figure(140)
clf
semilogy(0:length(result.erroriter)-1,result.erroriter,courbestyle{1},'LineWidth',2,'MarkerSize',10)
h=gca;
set(gca,'Fontsize',fsizetransp)
xlim([0,length(result.erroriter)])
xlabel('iteration','FontSize',fsizetransp)
ylabel('error','FontSize',fsizetransp)
saveas(gcf,'E:\REDACTION\matlab_results\Lshaperandomfield\Lshape_conv_pssub_dim8.eps','epsc2')


%%
GSD = GSDSOLVER('tol',1e-4,'pfixtol',1e-7,'pfixmax',6,'nbfoncmax',8,...
    'display',true,...
    'errorindicator','none','orthocrit',1e-10,'type','arnoldi',...
    'direct',true);

for i=1:8
    GSD = setparam(GSD,'nbfoncmax',i);
    qrada = solve(GSD,Ksto,fsto);
    erroragsd(i)=norm(qrada-qpc)/norm(qpc);
    errorAagsd(i)=norm(qrada-qpc,Ksto)/norm(qpc,Ksto);
end
qradasd = spectral_decomposition(qrada,'nbfoncmax',8,'pfixmax',7,'pfixtol',1e-10,'display');


%% arnolid avec 1 fonction de plus et selection
GSD = GSDSOLVER('tol',1e-4,'pfixtol',1e-7,'pfixmax',6,'nbfoncmax',8,...
    'display',true,'orthocrit',1e-11,...
    'errorindicator','none','type','arnoldi',...
    'direct',true,'inittype','one');

nbplus = 5
for i=1:8
    GSD = setparam(GSD,'nbfoncmax',i+nbplus);
    qrada = solve(GSD,Ksto,fsto);
    qradasd = spectral_decomposition(qrada,'nbfoncmax',i,'pfixmax',7,'pfixtol',1e-10);
    erroragsdtri(i)=norm(qradasd-qpc)/norm(qpc);
    errorAagsdtri(i)=norm(qradasd-qpc,Ksto)/norm(qpc,Ksto);
end

figure(1444)
clf
nbfonc=8;
semilogy(errorpsgsd,courbestyle{2},'LineWidth',2,'MarkerSize',10)
hold on
semilogy( erroragsdtri,courbestyle{3},'LineWidth',2,'MarkerSize',10)
h=gca;
set(gca,'Fontsize',fsizetransp)
xlabel('order M','FontSize',fsizetransp)
ylabel('error','FontSize',fsizetransp)
legend('SI-GSD',['A^{M+' num2str(nbplus) '}-GSD'])
xlim([1,nbfonc])
saveas(gcf,['E:\REDACTION\matlab_results\Lshaperandomfield\Lshape_error_agsdMplus' num2str(nbplus) '_psgsd.eps'],'epsc2')

figure(1445)
clf
nbfonc=8;
semilogy(errorApsgsd,courbestyle{2},'LineWidth',2,'MarkerSize',10)
hold on
semilogy( errorAagsdtri,courbestyle{3},'LineWidth',2,'MarkerSize',10)
h=gca;
set(gca,'Fontsize',fsizetransp)
xlabel('order M','FontSize',fsizetransp)
ylabel('error A ','FontSize',fsizetransp)
legend('SI-GSD',['A^{M+' num2str(nbplus) '}-GSD'])
xlim([1,nbfonc])
saveas(gcf,['E:\REDACTION\matlab_results\Lshaperandomfield\Lshape_errorA_agsdMplus' num2str(nbplus) '_psgsd.eps'],'epsc2')



%%
figure(14)
clf
nbfonc=8;
semilogy(errorqpcsd,courbestyle{1},'LineWidth',2,'MarkerSize',10)
hold on
semilogy( errorpsgsd,courbestyle{2},'LineWidth',2,'MarkerSize',10)
h=gca;
set(gca,'Fontsize',fsizetransp)
xlabel('order M','FontSize',fsizetransp)
ylabel('error','FontSize',fsizetransp)
legend('SD','SI-GSD')
xlim([1,nbfonc])
saveas(gcf,'E:\REDACTION\matlab_results\Lshaperandomfield\Lshape_error_sd_psgsd.eps','epsc2')

figure(144)
clf
nbfonc=8;
semilogy(errorpsgsd,courbestyle{2},'LineWidth',2,'MarkerSize',10)
hold on
semilogy( erroragsd,courbestyle{3},'LineWidth',2,'MarkerSize',10)
h=gca;
set(gca,'Fontsize',fsizetransp)
xlabel('order M','FontSize',fsizetransp)
ylabel('error','FontSize',fsizetransp)
legend('SI-GSD','A-GSD')
xlim([1,nbfonc])
saveas(gcf,'E:\REDACTION\matlab_results\Lshaperandomfield\Lshape_error_agsd_psgsd.eps','epsc2')


figure(155)
clf
nbfonc=8;
semilogy(errorApsgsd,courbestyle{2},'LineWidth',2,'MarkerSize',10)
hold on
semilogy( errorAagsd,courbestyle{3},'LineWidth',2,'MarkerSize',10)
h=gca;
set(gca,'Fontsize',fsizetransp)
xlabel('order M','FontSize',fsizetransp)
ylabel('error A','FontSize',fsizetransp)
legend('SI-GSD','A-GSD')
xlim([1,nbfonc])
saveas(gcf,'E:\REDACTION\matlab_results\Lshaperandomfield\Lshape_errorA_agsd_psgsd.eps','epsc2')



figure(15)
clf
nbfonc=8;
semilogy(errorAqpcsd,courbestyle{1},'LineWidth',2,'MarkerSize',10)
hold on
semilogy( errorApsgsd,courbestyle{2},'LineWidth',2,'MarkerSize',10)
h=gca;
set(gca,'Fontsize',fsizetransp)
xlabel('order M','FontSize',fsizetransp)
ylabel('error A','FontSize',fsizetransp)
legend('SD','SI-GSD')
xlim([1,nbfonc])
saveas(gcf,'E:\REDACTION\matlab_results\Lshaperandomfield\Lshape_errorA_sd_psgsd.eps','epsc2')


%%

figure(49)
clf
plotV(qpcsd,1:8,S,'nl',3,'manutext',{{1.2,0.3},'U','FontSize',18},'fact',.95)
%saveas(gcf,'E:\REDACTION\matlab_results\Lshaperandomfield\Lshape_sd_spacefun_sorted.jpg','jpeg')
figure(50)
clf
plotV(qradps,1:8,S,'nl',3,'manutext',{{1.2,0.3},'U','FontSize',18},'fact',.95)
%saveas(gcf,'E:\REDACTION\matlab_results\Lshaperandomfield\Lshape_psgsd_spacefun_sorted.jpg','jpeg')
figure(51)
clf
qradsd = spectral_decomposition(qrad,'nbfoncmax',10,'pfixmax',7,'pfixtol',1e-10,'display');
plotV(qradpssd,1:8,S,'nl',3,'manutext',{{1.2,0.3},'U','FontSize',18},'fact',.95)
%saveas(gcf,'E:\REDACTION\matlab_results\Lshaperandomfield\Lshape_psgsd_spacefun.jpg','jpeg')

figure(52)
clf
plotV(qrada,1:8,S,'nl',3,'manutext',{{1.2,0.3},'U','FontSize',18},'fact',.95)
%saveas(gcf,'E:\REDACTION\matlab_results\Lshaperandomfield\Lshape_agsd_spacefun.jpg','jpeg')
figure(53)
clf
qradasd = spectral_decomposition(qrada,'nbfoncmax',10,'pfixmax',7,'pfixtol',1e-10,'display');
plotV(qradasd,1:8,S,'nl',3,'manutext',{{1.2,0.3},'U','FontSize',18},'fact',.95)
%saveas(gcf,'E:\REDACTION\matlab_results\Lshaperandomfield\Lshape_agsd_spacefun_sorted.jpg','jpeg')




%%
%%
for i=1:getm(qpcsd)
    errorpsgsd(i)=norm(truncate(qradsd,1:i)-qpc)/norm(qpc);
    errorApsgsd(i)=norm(truncate(qradsd,1:i)-qpc,Ksto)/norm(qpc,Ksto);
end


%% Manufacture

qpcsd = spectral_decomposition(qpc,'nbfoncmax',4,'pfixmax',7,'pfixtol',1e-10);
qpcsd=truncate(qpcsd,[1,2]);

csto = Ksto*qpcsd;
clear errormanu_arnoldi
for i=1:5
    GSD = GSDSOLVER('tol',1e-4,'nbfoncmax',i+0,...
        'display',true,...
        'errorindicator','reference',...
        'type','arnoldi','direct',false,'orthocrit',1e-10);
    [qrada,result_arnoldi]=solve(GSD,Ksto,csto,[],'reference',qpcsd);
    qradasd = spectral_decomposition(qrada,'nbfoncmax',i,'pfixmax',7,'pfixtol',1e-10);
    errormanu_arnoldi(i)=norm(qradasd-qpcsd)/norm(qpcsd)
end

figure(30)
clf
plot(getV(qpcsd,1),S)
axis off
title('U_1','FontSize',24)
saveas(gcf,'E:\REDACTION\matlab_results\Lshaperandomfield\Lshape_manufacture_mode1_reference.jpg','jpeg')
saveas(gcf,'E:\REDACTION\matlab_results\Lshaperandomfield\Lshape_manufacture_mode1_reference.eps','epsc2')

figure(31)
clf
multipdfplot(getL(qpcsd),'nbs',1e5)
axis off
view(2)
%title('PDF (\lambda_1,\lambda_2)','FontSize',24)
saveas(gcf,'E:\REDACTION\matlab_results\Lshaperandomfield\Lshape_manufacture_pdflambda12_reference.jpg','jpeg')
saveas(gcf,'E:\REDACTION\matlab_results\Lshaperandomfield\Lshape_manufacture_pdflambda12_reference.eps','epsc2')


figure(32)
clf
plot(getV(qpcsd,2),S)
axis off
title('U_2','FontSize',24)
saveas(gcf,'E:\REDACTION\matlab_results\Lshaperandomfield\Lshape_manufacture_mode2_reference.jpg','jpeg')
saveas(gcf,'E:\REDACTION\matlab_results\Lshaperandomfield\Lshape_manufacture_mode2_reference.eps','epsc2')



clear errormanu_powersubspace
for i=0:4
    GSD = GSDSOLVER('tol',1e-4,'pfixtol',1e-17,'pfixmax',i,'nbfoncmax',2,...
        'display',true,...
        'errorindicator','residual','type','powersubspace','direct',false);
    [qrad,result_powersubspace]=solve(GSD,Ksto,csto,[],'reference',qpcsd);
    errormanu_powersubspace(i+1)=norm(qrad-qpcsd)/norm(qpcsd)
end

