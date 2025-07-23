%% Definition probleme
%  P =  70   126   210
dim=2
fsize = 14;
choixmail=2;
choixboun=2;

stdk = 0.1925;
stdsec = 0.1925;
% stdk=.3;
% stdsec=.3;


mk0 = 3;
mk1 = 1.5;
mf = 6;
mg = 2.25;
sk0 = stdk;
sk1 = stdk;
sf = stdsec;
sg = stdsec;

typemesh = 'fin'
psol = 3;

switch dim
    case 2
        r=15;
        P = POINT([ 0,0 ; 1,0 ; 0,2 ; 1,2 ; 1,1 ; 2,1 ; 2,2 ; 0,1 ]);
        if choixmail==2
            S1 = cast2matlab_model('PLANE',['C:\PROGRAMMES\CASTEM\EXAMPLES\Lshape_domain1_' typemesh '.txt']);
            S2 = cast2matlab_model('PLANE',['C:\PROGRAMMES\CASTEM\EXAMPLES\Lshape_domain2_' typemesh '.txt']);
            S = union(S1,S2);
            S=splitgroupelem(S,5000);
            if choixboun==1
                S = concatgroupelem(S);
            end
            S = createddlnode(S,DDL('u'),DDL('g'));
            if choixboun==1
                S = addcl(S,create_boundary(S),'u',0);
            else
                S = addcl(S,LINE(P(6),P(7)),'u',0);
            end
        else
            S1 = mesh(DOMAIN(2,P(1),P(5)),r,r);
            S2 = mesh(DOMAIN(2,P(8),P(4)),r,r);
            S3 = mesh(DOMAIN(2,P(5),P(7)),r,r);
            S = union(S1,S2);
            S = union(S,S3);
            S = concatgroupelem(S);
            S = convertelem(S,'TRI3')
            S = createddlnode(S,DDL('u'),DDL('g'));
            if choixboun==1
                S = addcl(S,create_boundary(S),'u',0);
            else
                S = addcl(S,LINE(P(6),P(7)),'u',0);
            end
        end
        
    case 1
        P = POINT([0;1]);
        S = mesh(DOMAIN(1,P(1),P(2)),50);
        S = createddlnode(S,DDL('u'),DDL('g'));
        S = addcl(S,P,'u',0);
        
end
%


choixsto=2;
if choixsto==1
    RV{1}=RVUNIFORM(0.7,1.3);
    clear X;
    PCM = PCMODEL(RV,'order',psol,'pcg');
    X{1} = PCM{1};
    
    X{2} = 1;
elseif choixsto==2
    jjj0=1;
    jjj1=1;
    nbsec=2;
    detrhs=0;
    RV = RANDVARS();
    for i=1:jjj0+jjj1+(nbsec*(detrhs~=1))
        RV{i} = RVUNIFORM(-sqrt(3),sqrt(3));
    end
    PCM = PCMODEL(RV,'order',psol,'pcg');
    PC =getPC(PCM);
    clear X;
    X{1}=mk0;
    if jjj0>0
        X{1}=mk0*(1+sk0*sum(PCM(1:jjj0),1)/jjj0);
    end
    X{2}=mk1;
    if jjj1>0
        X{2}=mk1*(1+sk1*sum(PCM(jjj0+(1:jjj1)),1)/jjj1);
    end
    if detrhs==1
        X{3}=mf;
        if nbsec==2
            X{4}=mg;
        end
    else
        X{3}=mf*(1+sf*full(PCM(jjj0+jjj1+1)));
        if nbsec==2
            X{4}=mg*(1+sg*full(PCM(jjj0+jjj1+2)));
        end
    end
    
end

thermique_nonlin_forms

%
ka=mean(X{1});
kn=mean(X{2});
if nbsec==2
    kl={mean(X{3}),mean(X{4})};
else
    kl=mean(X{3});
end
u = GSDthermiquesolvedet(S,1e-12,ka,kn,kl,'display');

if nbsec==1
    ulin = (ka*a{S}(:,:))\(kl*l{S}(:));
else
    ulin = (ka*a{S}(:,:))\(kl{1}*l{S}(:)+kl{2}*l2{S}(:));
end

norm(u-ulin)/norm(u)
figure(1)
clf
if dim==2
    subplot(1,3,1)
    plot(FENODEFIELD(ulin),S,'surface');
    title('lineaire')
    ax=axis;
    subplot(1,3,2)
    plot(FENODEFIELD(u),S,'surface');
    title('nonlineaire')
    axis(ax)
    subplot(1,3,3)
    plot(FENODEFIELD(-u+ulin),S,'surface');
    title('difference')
    axis(ax)
else
    plot(getcoord(getnode(S)),unfreevector(S,u),'r')
    hold on
    plot(getcoord(getnode(S)),unfreevector(S,ulin),'b')
    legend('non-lineaire','lineaire')
end

figure(3)
plot(u,S,'courbe','Color','r')
hold on

%% solution, de reference
cas=2
[qref,resultref] = thermique_solvereference(1e-3,S,PC,X,cas);
if cas==1
    fich=['diff_ref_typemesh_' typemesh '_stdk_' num2str(stdk) '_stdsec_' num2str(stdsec) '_mk0_' num2str(mk0) '_mk1_' num2str(mk1) '_mf_' num2str(mf) '_mg_' num2str(mg) '_p' num2str(psol) '.mat'];
else
    fich=['diff_ref2_typemesh_' typemesh '_stdk_' num2str(stdk) '_stdsec_' num2str(stdsec) '_mk0_' num2str(mk0) '_mk1_' num2str(mk1) '_mf_' num2str(mf) '_mg_' num2str(mg) '_p' num2str(psol) '.mat'];
end
% save(fich,'qref','resultref')
uMrefL2 = spectral_decomposition(qref,'nbfoncmax',12);
%% Monte-Carlo
Q=1000;
umc = zeros(getnbddlfree(S),Q);
for q=1:Q
    pourcentage(q,Q,100)
    for i=1:4
        Xq{i} = random(X{i});
    end
    umc(:,q) = GSDthermiquesolvedet(S,1e-8,Xq{1},Xq{2},{Xq{3},Xq{4}});
end

%%
%
nbfun = 12;
pfixiter = 3;
tolfinal = 1e-6;
epspfix = 2e-1;
toliter=epspfix/100;
update=0;
subspace=0;
ortho=0;

result=struct();
result.nbiter = zeros(1,nbfun);
result.nbitercumul = zeros(1,nbfun);

uM = PCRADIALMATRIX([size(u)],PC);
uM0=uM;
stoU = zeros(size(u,1),0);
stol = zeros(0,PC);
clear errGSD
clear errres

fun_G = @(W,L0) thermique_nonlin_G(W,tolfinal,S,PC,X,L0);
fun_GM = @(U,uM) thermique_nonlin_GM(U,uM,toliter,S,PC,X);
fun_F = @(L) thermique_nonlin_F(L,tolfinal,S,X);
fun_FM = @(l,uM) thermique_nonlin_FM(l,uM,toliter,S,X);
fun_res = @(uM) thermique_nonlin_res(uM,S,X,tolfinal/10);
fun_reserror = @(uM) norm(fun_res(uM))/norm(fun_res([]));

clock0=clock;
for M=1:nbfun
    fprintf('-----------------------------\n GSD couple #%d\n-----------------------------\n',M)
    lam0=allones(PC);
    alpha0 = norm(lam0);
    lam0 = lam0/alpha0;
    
    fprintf('-----------------------------\n POWER ITERATIONS \n-----------------------------\n')
    for iterpfix=1:pfixiter
        
        U = fun_FM(lam0,uM);
        
        U = U/norm(U);
        lam = fun_GM(U,uM);
        alpha = norm(lam);
        err = abs(alpha-alpha0)/(1/2*(alpha+alpha0));
        alpha0 = alpha;
        lam0=lam;
        
        fprintf('POWER iteration #%d : erreur %.3d\n',iterpfix,err);
        if err<epspfix && iterpfix>1
            break
        end
    end
    
    result.nbiter(M)=iterpfix;
    result.nbitercumul(M)=iterpfix;
    if M>1
        result.nbitercumul(M)=result.nbitercumul(M-1)+iterpfix;
        result.errL2uM(M) = norm(U.*lam)/norm(uM);
    else
        result.nbitercumul(M)=iterpfix;
        result.errL2uM(M) = 1;
    end
    
    stoU = [stoU,U];
    stol = [stol;lam];
    
    if ortho~=0
        [stoU,R] = qr(stoU,0);
        stol = R*stol;
    end
    
    if update==0 || M==1
        uMtemp = PCRADIALMATRIX(stoU,[size(u)],stol);
        deltauM = expand(uM-uMtemp);
        % deltauM = spectral_decomposition(uM-uMtemp,'tol',1e-15);
        errGSD(M) = norm(deltauM)/norm(uMtemp);
        uM = uMtemp;
    elseif update==1
        for kkk=0:subspace
            fprintf('Reactualisation des va ...\n')
            
            stol = fun_G(stoU,stol);
            
            if kkk<subspace
                fprintf('Reactualisation des vecteurs ...\n')
                stoU = fun_F(stol);
                [stoU,R] = qr(stoU,0);
                stol = R*stol;
            end
        end
        uMtemp = PCRADIALMATRIX(stoU,[size(u)],stol);
        % deltauM = expand(uM-uMtemp);
        % deltauM = spectral_decomposition(uM-uMtemp,'tol',1e-15);
        % errGSD(M) = norm(deltauM)/norm(uMtemp);
        uM = uMtemp;
    end
    % fprintf('\n GSD : number of functions = %d , erreur %.3d\n',M,errGSD(M));
    
    errres(getm(uM)) = fun_reserror(uM);
    fprintf(' GSD : number of functions = %d , erreur residu %.3d\n',getm(uM),errres(getm(uM)));
    result.time(getm(uM)) = etime(clock,clock0);
end

if update==1
    errresPUGSD = errres;
    result.errres = errres;
    resultPUGSD = result;
    uMPUGSD=uM;
    
    fich=['diff_PUGSD_typemesh_' typemesh '_stdk_' num2str(stdk) '_stdsec_' num2str(stdsec) '_mk0_' num2str(mk0) '_mk1_' num2str(mk1) '_mf_' num2str(mf) '_mg_' num2str(mg) '_p' num2str(psol) '.mat'];
    % save(fich,'uMPUGSD','resultPUGSD','errresPUGSD')
    
else
    errresPGSD = errres;
    result.errres = errres;
    resultPGSD = result;
    uMPGSD=uM;
    fich=['diff_PGSD_typemesh_' typemesh '_stdk_' num2str(stdk) '_stdsec_' num2str(stdsec) '_mk0_' num2str(mk0) '_mk1_' num2str(mk1) '_mf_' num2str(mf) '_mg_' num2str(mg) '_p' num2str(psol) '.mat'];
    % save(fich,'uMPGSD','resultPGSD','errresPGSD')
end
etime(clock,clock0)


upc  = unfreevector(S,uM);

%% ARNOLDI

nbfun = 12;
tolfinal = 1e-5;
update=0;
ortho=0;
orthocrit=1e-9;
result=struct();
testerror=1;
uM = PCRADIALMATRIX(size(u),PC);
uM0=uM;
clear errGSD
clear errres

fun_G = @(W,L0) thermique_nonlin_G(W,tolfinal,S,PC,X,L0);
fun_GM = @(U,uM) thermique_nonlin_GM(U,uM,tolfinal,S,PC,X);
fun_F = @(L) thermique_nonlin_F(L,tolfinal,S,X);
fun_FM = @(l,uM) thermique_nonlin_FM(l,uM,tolfinal,S,X);
fun_res = @(uM) thermique_nonlin_res(uM,S,X,tolfinal/10);
fun_reserror = @(uM) norm(fun_res(uM))/norm(fun_res([]));

fprintf('-----------------------------\n ARNOLDI ITERATIONS \n-----------------------------\n')

clock0=clock;
for rest = 0:0
    lam0=allones(PC);
    alpha0 = norm(lam0);
    lam0 = lam0/alpha0;
    stoU = zeros(size(u,1),0);
    stol = zeros(0,PC);
    
    fprintf('Restart %d \n--------\n',rest)
    
    for jj=1:nbfun-getm(uM)
        
        U = fun_FM(lam0,uM);
        U = U/norm(U);
        for ll=1:size(stoU,2)
            U = U - (U'*stoU(:,ll))*stoU(:,ll);
        end
        
        fprintf('\nresidu = %.3e',full(norm(U)))
        
        if norm(U)<orthocrit
            disp(' -> break')
            break
        end
        U = U /norm(U) ;
        stoU=[stoU,U];
        
        fprintf(' -> %d fonctions\n',size(stoU,2))
        
        lam = fun_GM(U,uM);
        lam0=lam;
        
        if testerror
            
            if size(stoU,2)==1
                uMtemp = PCRADIALMATRIX(U,size(u),lam);
            else
                stol = fun_G(stoU,zeros(size(stoU,2),PC));
                uMtemp = PCRADIALMATRIX(stoU,size(u),stol);
            end
            errres(getm(uMtemp)) = fun_reserror(uMtemp);
        end
        
    end
    
    fprintf('\nCalcul des va ...\n')
    stol = fun_G(stoU,zeros(size(stoU,2),PC));
    
    uM = PCRADIALMATRIX(stoU,size(u),stol);
    
    errres(getm(uM)) = fun_reserror(uM);
    fprintf('\nGSD : number of functions = %d , erreur residu %.3d\n',getm(uM),errres(getm(uM)));
    result.time(getm(uM)) = etime(clock,clock0);
    
end


errresAGSD = errres;
result.errres = errres;
resultAGSD = result;
uMAGSD=uM;

fich=['diff_AGSD_typemesh_' typemesh '_stdk_' num2str(stdk) '_stdsec_' num2str(stdsec) '_mk0_' num2str(mk0) '_mk1_' num2str(mk1) '_mf_' num2str(mf) '_mg_' num2str(mg) '_p' num2str(psol) '.mat'];
save(fich,'uMAGSD','resultAGSD','errresAGSD')
etime(clock,clock0)

upc  = unfreevector(S,uM);


%%

pp=3;
typemesh='fin';
tolaff = 5e-3
leg={};
figure(10)
clf
fich=['diff_ref2_typemesh_' typemesh '_stdk_' num2str(stdk) '_stdsec_' num2str(stdsec) '_mk0_' num2str(mk0) '_mk1_' num2str(mk1) '_mf_' num2str(mf) '_mg_' num2str(mg) '_p' num2str(pp) '.mat'];
load(fich)
rep=find(resultref.error>=tolaff);if length(rep)<length(resultref.error);rep=[rep,rep(end)+1];end
semilogy(resultref.time(rep),resultref.error(rep),['ks-'])
leg = [leg , {['ref.']}];
hold on
fich=['diff_PGSD_typemesh_' typemesh '_stdk_' num2str(stdk) '_stdsec_' num2str(stdsec) '_mk0_' num2str(mk0) '_mk1_' num2str(mk1) '_mf_' num2str(mf) '_mg_' num2str(mg) '_p' num2str(pp) '.mat'];
load(fich)
rep=find(resultPGSD.errres>=tolaff);if length(rep)<length(resultPGSD.errres);rep=[rep,rep(end)+1];end
semilogy(resultPGSD.time(rep),resultPGSD.errres(rep),['b*-'])
leg = [leg , {['algo 1']}];
fich=['diff_PUGSD_typemesh_' typemesh '_stdk_' num2str(stdk) '_stdsec_' num2str(stdsec) '_mk0_' num2str(mk0) '_mk1_' num2str(mk1) '_mf_' num2str(mf) '_mg_' num2str(mg) '_p' num2str(pp) '.mat'];
load(fich)
rep=find(resultPUGSD.errres>=tolaff);if length(rep)<length(resultPUGSD.errres);rep=[rep,rep(end)+1];end
semilogy(resultPUGSD.time(rep),resultPUGSD.errres(rep),['rv-'])
leg = [leg , {['algo 2']}];
set(gca,'FontSize',fsize)
legend(leg{:})
ylabel('Residual Error')
xlabel('Time (s)')

%%
scanpp = 3:5;
ppline = {'-','--','-.','.'};
% % % %
tolaff = 5e-3
typemesh = 'hyper_fin';
figure(10)
clf
leg = {};
for pp=scanpp
    try
        fich=['diff_ref2_typemesh_' typemesh '_stdk_' num2str(stdk) '_stdsec_' num2str(stdsec) '_mk0_' num2str(mk0) '_mk1_' num2str(mk1) '_mf_' num2str(mf) '_mg_' num2str(mg) '_p' num2str(pp) '.mat'];
        load(fich)
        rep=find(resultref.error>=tolaff);if length(rep)<length(resultref.error);rep=[rep,rep(end)+1];end
        semilogy(resultref.time(rep),resultref.error(rep),['ks' ppline{pp==scanpp}])
        leg = [leg , {['ref., N_0=' num2str(pp)]}];
        hold on
    catch
        warning(fich)
    end
    
    try
        fich=['diff_PGSD_typemesh_' typemesh '_stdk_' num2str(stdk) '_stdsec_' num2str(stdsec) '_mk0_' num2str(mk0) '_mk1_' num2str(mk1) '_mf_' num2str(mf) '_mg_' num2str(mg) '_p' num2str(pp) '.mat'];
        load(fich)
        rep=find(resultPGSD.errres>=tolaff);if length(rep)<length(resultPGSD.errres);rep=[rep,rep(end)+1];end
        semilogy(resultPGSD.time(rep),resultPGSD.errres(rep),['b*' ppline{pp==scanpp} ])
        leg = [leg , {['algo 1, N_0=' num2str(pp)]}];
        hold on
    catch
        warning(fich)
    end
    try
        fich=['diff_PUGSD_typemesh_' typemesh '_stdk_' num2str(stdk) '_stdsec_' num2str(stdsec) '_mk0_' num2str(mk0) '_mk1_' num2str(mk1) '_mf_' num2str(mf) '_mg_' num2str(mg) '_p' num2str(pp) '.mat'];
        load(fich)
        rep=find(resultPUGSD.errres>=tolaff);if length(rep)<length(resultPUGSD.errres);rep=[rep,rep(end)+1];end
        semilogy(resultPUGSD.time(rep),resultPUGSD.errres(rep),['rv' ppline{pp==scanpp} ])
        leg = [leg , {['algo 2, N_0=' num2str(pp)]}];
        hold on
    catch
        warning(fich)
    end
end
set(gca,'FontSize',fsize)
legend(leg{:})
ylabel('Residual Error')
xlabel('Time (s)')

%%
scantypemesh = {'gros','fin','super_fin','hyper_fin'};
scann=[178,368,726,1431]
scantypemesh=scantypemesh(2:4)
scann=scann(2:4)
ppline = {'-','--','-.',':'};
pp=3;
% % % %
tolaff = 5e-3
figure(10)
clf
leg = {};
for kk=1:length(scantypemesh)
    typemesh=scantypemesh{kk};
    try
        fich=['diff_ref2_typemesh_' typemesh '_stdk_' num2str(stdk) '_stdsec_' num2str(stdsec) '_mk0_' num2str(mk0) '_mk1_' num2str(mk1) '_mf_' num2str(mf) '_mg_' num2str(mg) '_p' num2str(pp) '.mat'];
        load(fich)
        rep=find(resultref.error>=tolaff);if length(rep)<length(resultref.error);rep=[rep,rep(end)+1];end
        semilogy(resultref.time(rep),resultref.error(rep),['ks' ppline{kk}])
        leg = [leg , {['ref., N_x=' num2str(scann(kk))]}];
        hold on
    catch
        warning(fich)
    end
    
    try
        fich=['diff_PGSD_typemesh_' typemesh '_stdk_' num2str(stdk) '_stdsec_' num2str(stdsec) '_mk0_' num2str(mk0) '_mk1_' num2str(mk1) '_mf_' num2str(mf) '_mg_' num2str(mg) '_p' num2str(pp) '.mat'];
        load(fich)
        rep=find(resultPGSD.errres>=tolaff);if length(rep)<length(resultPGSD.errres);rep=[rep,rep(end)+1];end
        semilogy(resultPGSD.time(rep),resultPGSD.errres(rep),['b*' ppline{kk} ])
        leg = [leg , {['algo 1, N_x=' num2str(scann(kk))]}];
        hold on
    catch
        
        warning(fich)
    end
    try
        fich=['diff_PUGSD_typemesh_' typemesh '_stdk_' num2str(stdk) '_stdsec_' num2str(stdsec) '_mk0_' num2str(mk0) '_mk1_' num2str(mk1) '_mf_' num2str(mf) '_mg_' num2str(mg) '_p' num2str(pp) '.mat'];
        load(fich)
        rep=find(resultPUGSD.errres>=tolaff);if length(rep)<length(resultPUGSD.errres);rep=[rep,rep(end)+1];end
        semilogy(resultPUGSD.time(rep),resultPUGSD.errres(rep),['rv' ppline{kk} ])
        leg = [leg , {['algo 2, N_x=' num2str(scann(kk))]}];
        hold on
    catch
        warning(fich)
    end
end
set(gca,'FontSize',fsize)
legend(leg{:})
ylabel('Residual Error')
xlabel('Time (s)')


%%
scantypemesh = {'gros','fin','super_fin','hyper_fin'};
scann=[178,368,726,1431];
% scantypemesh = scantypemesh(2:4);
% scann=scann(2:4);
scanp=[3:6];
scanP = factorial(4+scanp)/factorial(4)./factorial(scanp);
tol = 5e-2;
timeref = zeros(length(scann),length(scanp));
timeref2 = zeros(length(scann),length(scanp));
timepgsd = timeref;
timepugsd = timeref;
Mpgsd = timeref;
Mpugsd = timeref;
errpgsd = timeref;
errpugsd = timeref;

for kk=1:length(scantypemesh)
    typemesh=scantypemesh{kk};
    for ll=1:length(scanp)
        pp = scanp(ll);
        try
            fich=['diff_ref_typemesh_' typemesh '_stdk_' num2str(stdk) '_stdsec_' num2str(stdsec) '_mk0_' num2str(mk0) '_mk1_' num2str(mk1) '_mf_' num2str(mf) '_mg_' num2str(mg) '_p' num2str(pp) '.mat'];
            load(fich)
            rep=find(resultref.error<=tol);rep=rep(1);
            timeref(kk,ll) = resultref.time(rep);
        catch
            warning(fich)
        end
        try
            fich=['diff_ref2_typemesh_' typemesh '_stdk_' num2str(stdk) '_stdsec_' num2str(stdsec) '_mk0_' num2str(mk0) '_mk1_' num2str(mk1) '_mf_' num2str(mf) '_mg_' num2str(mg) '_p' num2str(pp) '.mat'];
            load(fich)
            rep=find(resultref.error<=tol);rep=rep(1);
            timeref2(kk,ll) = resultref.time(rep);
        catch
            warning(fich)
        end
        try
            fich=['diff_PGSD_typemesh_' typemesh '_stdk_' num2str(stdk) '_stdsec_' num2str(stdsec) '_mk0_' num2str(mk0) '_mk1_' num2str(mk1) '_mf_' num2str(mf) '_mg_' num2str(mg) '_p' num2str(pp) '.mat'];
            load(fich)
            rep=find(resultPGSD.errres<=tol);rep=rep(1);
            timepgsd(kk,ll) = resultPGSD.time(rep);
            Mpgsd(kk,ll) = rep;
            errpgsd(kk,ll) = resultPGSD.errres(rep);
        catch
            warning(fich)
        end
        try
            fich=['diff_PUGSD_typemesh_' typemesh '_stdk_' num2str(stdk) '_stdsec_' num2str(stdsec) '_mk0_' num2str(mk0) '_mk1_' num2str(mk1) '_mf_' num2str(mf) '_mg_' num2str(mg) '_p' num2str(pp) '.mat'];
            load(fich)
            rep=find(resultPUGSD.errres<=tol);rep=rep(1);
            timepugsd(kk,ll) = resultPUGSD.time(rep);
            Mpugsd(kk,ll) = rep;
            errpugsd(kk,ll) = resultPUGSD.errres(rep);
            
        catch
            warning(fich)
        end
    end
end

% timeref(:,end)=[11,30,43,61]*10^3
if tol==1e-2
    timeref2(end-2:end,end)=[5,6.5,8.5]*10^3*1.1
elseif tol==2e-2
    timeref2(end-2:end,end)=[5,6.5,12]*10^3*0.8
elseif tol==5e-2
    timeref2(end-2:end,end)=[5,6.5,8.5]*10^3*0.6
else
    error()
end


timeref=timeref2;
scannP = scann(:)*scanP;
figure(16)
clf
figure(160)
clf
leg={};leg2={};
for kk=1:length(scanp)
    figure(16)
    loglog(scannP(:,kk),timeref(:,kk)./timepgsd(:,kk),getcourbestyles(kk,'Marker'))
    leg = [leg , {['N_0=' num2str(scanp(kk))]}];
    hold on
    figure(160)
    loglog(scannP(:,kk),timeref(:,kk)./timepugsd(:,kk),getcourbestyles(kk,'Marker'))
    hold on
    leg2 = [leg2 , {['N_0=' num2str(scanp(kk))]}];
end
figure(16)
set(gca,'FontSize',fsize)
legend(leg{:})
ylabel('Time Gain')
xlabel('N_x\times P')
figure(160)
set(gca,'FontSize',fsize)
legend(leg2{:})
ylabel('Time Gain')
xlabel('N_x\times P')


figure(17)
clf
figure(170)
clf
figure(18)
clf
leg={};
legref={};
figure(19)
clf
leg={};
for kk=1:length(scann)
    figure(17)
    loglog(scannP(kk,:),timeref(kk,:)./timepgsd(kk,:),getcourbestyles(kk,'Marker'))
    hold on
    figure(170)
    loglog(scannP(kk,:),timeref(kk,:)./timepugsd(kk,:),getcourbestyles(kk,'Marker'))
    hold on
    figure(19)
    loglog(scannP(kk,:),timepgsd(kk,:),getcourbestyles(kk,'Marker'))
    hold on
    figure(190)
    loglog(scannP(kk,:),timepugsd(kk,:),getcourbestyles(kk,'Marker'))
    hold on
    figure(18)
    loglog(scannP(kk,:),timeref(kk,:),getcourbestyles(kk,'Marker'))
    hold on
    leg = [leg , {['N_x=' num2str(scann(kk))]}];
    %leg = [leg , {['N_x=' num2str(scann(kk))]}];
    legref = [legref , {['N_x=' num2str(scann(kk))]}];
end
figure(17)
set(gca,'FontSize',fsize)
legend(leg{:})
ylabel('Time Gain')
xlabel('N_x\times P')
figure(170)
set(gca,'FontSize',fsize)
legend(leg{:})
ylabel('Time Gain')
xlabel('N_x\times P')
figure(18)
set(gca,'FontSize',fsize)
legend(legref{:})
ylabel('Time (s)')
xlabel('N_x\times P')
figure(19)
set(gca,'FontSize',fsize)
legend(leg{:})
ylabel('Time (s)')
xlabel('N_x\times P')
figure(190)
set(gca,'FontSize',fsize)
legend(leg{:})
ylabel('Time (s)')
xlabel('N_x\times P')


timeref
Mpgsd
errpgsd
timepgsd
Mpugsd
errpugsd
timepugsd
timegainpgsd = timeref./timepgsd
timegainpugsd = timeref./timepugsd

figure(18)
%%

figure(100)
clf
loglog(scannP,timeref)


%%

for kk=1:length(scantypemesh)
    typemesh=scantypemesh{kk};
    for ll=1:length(scanp)
        pp = scanp(ll);
        try
            fich=['diff_ref_typemesh_' typemesh '_stdk_' num2str(stdk) '_stdsec_' num2str(stdsec) '_mk0_' num2str(mk0) '_mk1_' num2str(mk1) '_mf_' num2str(mf) '_mg_' num2str(mg) '_p' num2str(pp) '.mat'];
            load(fich)
        end
        
        try
            fich=['diff_PGSD_typemesh_' typemesh '_stdk_' num2str(stdk) '_stdsec_' num2str(stdsec) '_mk0_' num2str(mk0) '_mk1_' num2str(mk1) '_mf_' num2str(mf) '_mg_' num2str(mg) '_p' num2str(pp) '.mat'];
            load(fich)
            for i=1:length(uMPGSD)
                errL2pgsd(i) = norm(truncate(uMPGSD,1:i)-uMPGSD)/norm(uMPGSD);
            end
        end
        
        rep=find(resultPGSD.errres<=tol);fprintf('(res) %s, p=%d, M=%d\n',typemesh,pp,rep(1))
        rep=find(errL2pgsd<=tol);fprintf('(L2 ) %s, p=%d, M=%d\n',typemesh,pp,rep(1))
        
    end
end



%% CONVERGENCE DE PUGSD ET PGSD
nbmodes=nbfun;
figure(44)
clf
semilogy(1:nbmodes,resultPGSD.errL2uM(1:nbmodes),'b--s','LineWidth',2)
hold on
semilogy(1:nbmodes,errresPGSD(1:nbmodes),'r-o','LineWidth',2)
set(gca,'FontSize',fsize)
set(gca,'Xtick',[1:nbmodes])
xlim([.8,nbmodes])
xlabel('M')
ylabel('|| R_M ||')
legend('algorithm 1, ||\lambda_MU_M||/||u_M||','algorithm 1, ||R_M||')
%%
nbmodes=nbfun;
figure(44)
clf
semilogy(1:nbmodes,resultPUGSD.errL2uM(1:nbmodes),'b--s','LineWidth',2)
hold on
semilogy(1:nbmodes,errresPUGSD(1:nbmodes),'r-o','LineWidth',2)
set(gca,'FontSize',fsize)
set(gca,'Xtick',[1:nbmodes])
xlim([.8,nbmodes])
xlabel('M')
ylabel('|| R_M ||')
legend('algorithm 1, ||\lambda_MU_M||/||u_M||','algorithm 1, ||R_M||')
%%
nbmodes=nbfun;
figure(44)
clf
semilogy(1:nbmodes,resultPGSD.errL2uM(1:nbmodes),'b-s','LineWidth',2)
hold on
semilogy(1:nbmodes,resultPUGSD.errL2uM(1:nbmodes),'b-->','LineWidth',2)
semilogy(1:nbmodes,errresPGSD(1:nbmodes),'r-o','LineWidth',2)
semilogy(1:nbmodes,errresPUGSD(1:nbmodes),'r--*','LineWidth',2)
set(gca,'FontSize',fsize)
set(gca,'Xtick',[1:nbmodes])
xlim([.8,nbmodes])
xlabel('M')
ylabel('|| R_M ||')
legend('algorithm 1, ||\lambda_MU_M||/||u_M||','algorithm 2, ||\lambda_MU_M||/||u_M||','algorithm 1, ||R_M||','algorithm 2, ||R_M||')
nbmodes=nbfun;

%% CONVERGENCE PGSD PUGSD AGSD
nbmodes=nbfun;
figure(44)
clf
semilogy(1:nbmodes,errresPGSD(1:nbmodes),'b-o','LineWidth',2)
hold on
semilogy(1:nbmodes,errresPUGSD(1:nbmodes),'r-.*','LineWidth',2)
semilogy(1:nbmodes,errresAGSD(1:nbmodes),'k-->','LineWidth',2)
set(gca,'FontSize',fsize)
set(gca,'Xtick',[1:nbmodes])
xlim([.8,nbmodes])
xlabel('M')
ylabel('|| R_M ||')
legend('P-GSD','PU-GSD','A-GSD')

%% CONVERGENCE PGSD PUGSD AGSD
nbmodes=nbfun;
figure(45)
clf
semilogy(resultPGSD.nbitercumul,errresPGSD(1:nbmodes),'b-o','LineWidth',2)
hold on
semilogy(resultPUGSD.nbitercumul,errresPUGSD(1:nbmodes),'r-.*','LineWidth',2)
semilogy(1:nbmodes,errresAGSD(1:nbmodes),'k-->','LineWidth',2)
set(gca,'FontSize',fsize)
% set(gca,'Xtick',[1:nbmodes])
% xlim([.8,nbmodes])
xlabel('Number of Nonlinear Deterministic Problems')
ylabel('|| R_M ||')
legend('P-GSD','PU-GSD','A-GSD')

%%
% load('diffusion_nonlin_PGSD_p4_M16.mat')
% load('diffusion_nonlin_PUGSD_p4_M16.mat')
load('diff_PGSD_stdk_0.1925_stdsec_0.1925_mk0_3_mk1_1.5_mf_6_mg_2.25_p5.mat')
load('diff_PUGSD_stdk_0.1925_stdsec_0.1925_mk0_3_mk1_1.5_mf_6_mg_2.25_p5.mat')
nbmodes=8;
for i=1:nbmodes
    resultPGSD.errL2(i)=norm(truncate(uMPGSD,i+1:getm(uMPGSD)))/norm(uMPGSD);
    resultPUGSD.errL2(i)=norm(truncate(uMPUGSD,i+1:getm(uMPUGSD)))/norm(uMPUGSD);
end
%%

figure(44)
clf
semilogy(1:nbmodes,resultPGSD.errL2(1:nbmodes),'b-s','LineWidth',2)
hold on
semilogy(1:nbmodes,resultPUGSD.errL2(1:nbmodes),'b-->','LineWidth',2)
semilogy(1:nbmodes,errresPGSD(1:nbmodes),'r-o','LineWidth',2)
semilogy(1:nbmodes,errresPUGSD(1:nbmodes),'r--*','LineWidth',2)
set(gca,'FontSize',fsize)
set(gca,'Xtick',[1:nbmodes])
set(gca,'Ytick',10.^[-5:1])
xlim([.8,nbmodes])
y=ylim;
ylim([y(1),6])

xlabel('M')
ylabel('|| R_M ||')
legend('algorithm 1, ||u_M -u_{ref}||/||u_{ref}||','algorithm 2, ||u_M -u_{ref}||/||u_{ref}||','algorithm 1, ||R_M||','algorithm 2, ||R_M||')
%% CONVERGENCE DE PUGSD ET PGSD
nbmodes=nbfun;
figure(44)
clf
semilogy(1:nbmodes,errresPGSD(1:nbmodes),'b--s','LineWidth',2)
hold on
semilogy(1:nbmodes,errresPUGSD(1:nbmodes),'r-o','LineWidth',2)
set(gca,'FontSize',fsize)
set(gca,'Xtick',[1:nbmodes])
xlim([.8,nbmodes])
xlabel('M')
ylabel('|| R_M ||')
legend('algorithm 1','algorithm 2')
%% INFLUENCE STD
msize=8;
nbmodes = length(errresPGSD);
figure(1000)
clf
load('diff_PGSD_stdk_0.1_stdsec_0.1_mk0_3_mk1_1.5_mf_6_mg_2.25_p4.mat')
semilogy(1:length(errresPGSD),errresPGSD,'ks-','LineWidth',2,'MarkerSize',msize)
hold on
load('diff_PGSD_stdk_0.1925_stdsec_0.1925_mk0_3_mk1_1.5_mf_6_mg_2.25_p4.mat')
semilogy(1:length(errresPGSD),errresPGSD,'bo-','LineWidth',2,'MarkerSize',msize)
load('diff_PGSD_stdk_0.3_stdsec_0.3_mk0_3_mk1_1.5_mf_6_mg_2.25_p4.mat')
semilogy(1:length(errresPGSD),errresPGSD,'r^-','LineWidth',2,'MarkerSize',msize)
legend('cov = 10%','cov = 20%','cov = 30%')
set(gca,'FontSize',12)
set(gca,'Xtick',[1:nbmodes])
xlim([.8,nbmodes])
xlabel('M')
ylabel('|| R_M ||')



%%
msize=8;
nbmodes = length(errresPUGSD);
figure(1000)
clf
load('diff_PUGSD_stdk_0.1_stdsec_0.1_mk0_3_mk1_1.5_mf_6_mg_2.25_p4.mat')
semilogy(1:length(errresPUGSD),errresPUGSD,'ks-','LineWidth',2,'MarkerSize',msize)
hold on
load('diff_PUGSD_stdk_0.1925_stdsec_0.1925_mk0_3_mk1_1.5_mf_6_mg_2.25_p4.mat')
semilogy(1:length(errresPUGSD),errresPUGSD,'bo-','LineWidth',2,'MarkerSize',msize)
load('diff_PUGSD_stdk_0.3_stdsec_0.3_mk0_3_mk1_1.5_mf_6_mg_2.25_p4.mat')
semilogy(1:length(errresPUGSD),errresPUGSD,'r^-','LineWidth',2,'MarkerSize',msize)
legend('cov = 10%','cov = 20%','cov = 30%')
set(gca,'FontSize',12)
set(gca,'Xtick',[1:nbmodes])
xlim([.8,nbmodes])
xlabel('M')
ylabel('|| R_M ||')
%% INFLUENCE P

msize=8;
nbmodes = length(errresPGSD);
figure(1000)
clf
load('diff_PGSD_stdk_0.1925_stdsec_0.1925_mk0_3_mk1_1.5_mf_6_mg_2.25_p4.mat')
semilogy(1:length(errresPGSD),errresPGSD,'ks-','LineWidth',2,'MarkerSize',msize)
hold on
load('diff_PGSD_stdk_0.1925_stdsec_0.1925_mk0_3_mk1_1.5_mf_6_mg_2.25_p5.mat')
semilogy(1:length(errresPGSD),errresPGSD,'bo-','LineWidth',2,'MarkerSize',msize)
load('diff_PGSD_stdk_0.1925_stdsec_0.1925_mk0_3_mk1_1.5_mf_6_mg_2.25_p6.mat')
semilogy(1:length(errresPGSD),errresPGSD,'r^-','LineWidth',2,'MarkerSize',msize)
legend('N_o=4','N_o=5','N_o=6')
set(gca,'FontSize',12)
set(gca,'Xtick',[1:nbmodes])
xlim([.8,nbmodes])
xlabel('M')
ylabel('|| R_M ||')

%%
msize=8;
nbmodes = length(errresPUGSD);
figure(1000)
clf
load('diff_PUGSD_stdk_0.1925_stdsec_0.1925_mk0_3_mk1_1.5_mf_6_mg_2.25_p4.mat')
semilogy(1:length(errresPUGSD),errresPUGSD,'ks-','LineWidth',2,'MarkerSize',msize)
hold on
load('diff_PUGSD_stdk_0.1925_stdsec_0.1925_mk0_3_mk1_1.5_mf_6_mg_2.25_p5.mat')
semilogy(1:length(errresPUGSD),errresPUGSD,'bo-','LineWidth',2,'MarkerSize',msize)
load('diff_PUGSD_stdk_0.1925_stdsec_0.1925_mk0_3_mk1_1.5_mf_6_mg_2.25_p6.mat')
semilogy(1:length(errresPUGSD),errresPUGSD,'r^-','LineWidth',2,'MarkerSize',msize)
legend('N_o=4','N_o=5','N_o=6')
set(gca,'FontSize',12)
set(gca,'Xtick',[1:nbmodes])
xlim([.8,nbmodes])
xlabel('M')
ylabel('|| R_M ||')


%%
msize=8;
nbmodes = length(errresPUGSD);
figure(1000)
clf
load('diff_PGSD_stdk_0.1_stdsec_0.1_mk0_3_mk1_1.5_mf_6_mg_2.25_p4.mat')
semilogy(1:length(errresPGSD),errresPGSD,'ks-','LineWidth',2,'MarkerSize',msize)
hold on
load('diff_PUGSD_stdk_0.1_stdsec_0.1_mk0_3_mk1_1.5_mf_6_mg_2.25_p4.mat')
semilogy(1:length(errresPUGSD),errresPUGSD,'ks--','LineWidth',2,'MarkerSize',msize)
hold on
load('diff_PGSD_stdk_0.1925_stdsec_0.1925_mk0_3_mk1_1.5_mf_6_mg_2.25_p4.mat')
semilogy(1:length(errresPGSD),errresPGSD,'bo-','LineWidth',2,'MarkerSize',msize)
load('diff_PUGSD_stdk_0.1925_stdsec_0.1925_mk0_3_mk1_1.5_mf_6_mg_2.25_p4.mat')
semilogy(1:length(errresPUGSD),errresPUGSD,'bo--','LineWidth',2,'MarkerSize',msize)
load('diff_PGSD_stdk_0.3_stdsec_0.3_mk0_3_mk1_1.5_mf_6_mg_2.25_p4.mat')
semilogy(1:length(errresPGSD),errresPGSD,'r^-','LineWidth',2,'MarkerSize',msize)
load('diff_PUGSD_stdk_0.3_stdsec_0.3_mk0_3_mk1_1.5_mf_6_mg_2.25_p4.mat')
semilogy(1:length(errresPUGSD),errresPUGSD,'r^--','LineWidth',2,'MarkerSize',msize)
legend('algo 1 , cov = 10%','algo 2 , cov = 10%','algo 1 , cov = 20%','algo 2 , cov = 20%','algo 1 , cov = 30%','algo 2 , cov = 30%')
set(gca,'FontSize',12)
set(gca,'Xtick',[1:nbmodes])
xlim([.8,nbmodes])
xlabel('M')
ylabel('|| R_M ||')

%% INFLUENCE NONLIN (MOYENNE K1)
msize=8;
nbmodes = length(errresPUGSD);
figure(1000)
clf
load('diff_PUGSD_stdk_0.1925_stdsec_0.1925_mk0_3_mk1_1.5_mf_6_mg_2.25_p4.mat')
semilogy(1:length(errresPUGSD),errresPUGSD,'bo-','LineWidth',2,'MarkerSize',msize)
hold on
load('diff_PUGSD_stdk_0.1925_stdsec_0.1925_mk0_3_mk1_0.5_mf_6_mg_2.25_p4.mat')
semilogy(1:length(errresPUGSD),errresPUGSD,'r>-','LineWidth',2,'MarkerSize',msize)
load('diff_PUGSD_stdk_0.1925_stdsec_0.1925_mk0_3_mk1_0.1_mf_6_mg_2.25_p4.mat')
semilogy(1:length(errresPUGSD),errresPUGSD,'ks-','LineWidth',2,'MarkerSize',msize)
load('diff_PUGSD_stdk_0.1925_stdsec_0.1925_mk0_3_mk1_0.01_mf_6_mg_2.25_p4.mat')
semilogy(1:length(errresPUGSD),errresPUGSD,'m^-','LineWidth',2,'MarkerSize',msize)
load('diff_PUGSD_stdk_0.1925_stdsec_0.1925_mk0_3_mk1_2.2204e-016_mf_6_mg_2.25_p4.mat')
semilogy(1:length(errresPUGSD),errresPUGSD,'k*-','LineWidth',2,'MarkerSize',msize)
legend('E(\kappa_1) = 1.5','E(\kappa_1) = 0.5','E(\kappa_1) = 0.1','E(\kappa_1) = 0.01','\kappa_1 = 0')
set(gca,'FontSize',12)
set(gca,'Xtick',[1:nbmodes])
xlim([.8,nbmodes])
xlabel('M')
ylabel('|| R_M ||')

%%
msize=8;
nbmodes = length(errresPGSD);
figure(1000)
clf
load('diff_PGSD_stdk_0.1925_stdsec_0.1925_mk0_3_mk1_1.5_mf_6_mg_2.25_p4.mat')
semilogy(1:length(errresPGSD),errresPGSD,'bo-','LineWidth',2,'MarkerSize',msize)
hold on
load('diff_PGSD_stdk_0.1925_stdsec_0.1925_mk0_3_mk1_0.5_mf_6_mg_2.25_p4.mat')
semilogy(1:length(errresPGSD),errresPGSD,'r>-','LineWidth',2,'MarkerSize',msize)
load('diff_PGSD_stdk_0.1925_stdsec_0.1925_mk0_3_mk1_0.1_mf_6_mg_2.25_p4.mat')
semilogy(1:length(errresPGSD),errresPGSD,'ks-','LineWidth',2,'MarkerSize',msize)
load('diff_PGSD_stdk_0.1925_stdsec_0.1925_mk0_3_mk1_0.01_mf_6_mg_2.25_p4.mat')
semilogy(1:length(errresPGSD),errresPGSD,'m^-','LineWidth',2,'MarkerSize',msize)
load('diff_PGSD_stdk_0.1925_stdsec_0.1925_mk0_3_mk1_2.2204e-016_mf_6_mg_2.25_p4.mat')
semilogy(1:length(errresPGSD),errresPGSD,'k*-','LineWidth',2,'MarkerSize',msize)
legend('E(\kappa_1) = 1.5','E(\kappa_1) = 0.5','E(\kappa_1) = 0.1','E(\kappa_1) = 0.01','\kappa_1 = 0')
set(gca,'FontSize',12)
set(gca,'Xtick',[1:nbmodes])
xlim([.8,nbmodes])
xlabel('M')
ylabel('|| R_M ||')


%%
msize=8;
figure(1000)
clf
load('diffusion_nonlin_PGSD_p4_M16.mat')
semilogy(resultPGSD.nbitercumul,errresPGSD,'b--s','LineWidth',2,'MarkerSize',msize)
load('diffusion_nonlin_PUGSD_p4_M16.mat')
hold on
semilogy(resultPUGSD.nbitercumul,errresPUGSD,'r-o','LineWidth',2,'MarkerSize',msize)
set(gca,'FontSize',fsize)
% set(gca,'Xtick',max(max()))
xlim([-3,max(max(resultPGSD.nbitercumul),max(resultPUGSD.nbitercumul))]+3)
xlabel('iterations')
ylabel('|| R_M ||')
legend('algorithm 1','algorithm 2')


%% ROBUSTESSE EN FONCTION DE LA PRECISION DES ITERATIONS PUISSANCES
msize=8;
figure(1000)
clf
load('diff_PGSD_eps_0.5p4.mat')
semilogy(resultPGSD.nbitercumul,errresPGSD,'ks-','LineWidth',2,'MarkerSize',msize)
hold on
load('diff_PGSD_eps_0.1p4.mat')
semilogy(resultPGSD.nbitercumul,errresPGSD,'bo-','LineWidth',2,'MarkerSize',msize)
load('diff_PGSD_eps_0.01p4.mat')
semilogy(resultPGSD.nbitercumul,errresPGSD,'r^-','LineWidth',2,'MarkerSize',msize)
load('diff_PGSD_eps_0.001p4.mat')
semilogy(resultPGSD.nbitercumul,errresPGSD,'k>-','LineWidth',2,'MarkerSize',msize)
legend('\epsilon_s = 5.10^{-1}','\epsilon_s = 10^{-1}','\epsilon_s = 10^{-2}','\epsilon_s = 10^{-3}')
set(gca,'FontSize',12)
xlabel('iterations')
ylabel('|| R_M ||')
figure(1001)
clf
load('diff_PGSD_eps_0.5p4.mat')
semilogy(1:length(errresPGSD),errresPGSD,'ks-','LineWidth',2,'MarkerSize',msize)
hold on
load('diff_PGSD_eps_0.1p4.mat')
semilogy(1:length(errresPGSD),errresPGSD,'bo-','LineWidth',2,'MarkerSize',msize)
load('diff_PGSD_eps_0.01p4.mat')
semilogy(1:length(errresPGSD),errresPGSD,'r^-','LineWidth',2,'MarkerSize',msize)
load('diff_PGSD_eps_0.001p4.mat')
semilogy(1:length(errresPGSD),errresPGSD,'k>-','LineWidth',2,'MarkerSize',msize)
legend('\epsilon_s = 5.10^{-1}','\epsilon_s = 10^{-1}','\epsilon_s = 10^{-2}','\epsilon_s = 10^{-3}')
set(gca,'FontSize',12)
xlabel('M')
ylabel('|| R_M ||')

%%



msize=8;
figure(1000)
clf
load('diff_PUGSD_eps_0.5p4.mat')
semilogy(resultPUGSD.nbitercumul,errresPUGSD,'ks-','LineWidth',2,'MarkerSize',msize)
hold on
load('diff_PUGSD_eps_0.1p4.mat')
semilogy(resultPUGSD.nbitercumul,errresPUGSD,'bo-','LineWidth',2,'MarkerSize',msize)
load('diff_PUGSD_eps_0.01p4.mat')
semilogy(resultPUGSD.nbitercumul,errresPUGSD,'r^-','LineWidth',2,'MarkerSize',msize)
load('diff_PUGSD_eps_0.001p4.mat')
semilogy(resultPUGSD.nbitercumul,errresPUGSD,'k>-','LineWidth',2,'MarkerSize',msize)
legend('\epsilon_s = 5.10^{-1}','\epsilon_s = 10^{-1}','\epsilon_s = 10^{-2}','\epsilon_s = 10^{-3}')
set(gca,'FontSize',12)
xlabel('iterations')
ylabel('|| R_M ||')
figure(1001)
clf
load('diff_PUGSD_eps_0.5p4.mat')
semilogy(1:length(errresPUGSD),errresPUGSD,'ks-','LineWidth',2,'MarkerSize',msize)
hold on
load('diff_PUGSD_eps_0.1p4.mat')
semilogy(1:length(errresPUGSD),errresPUGSD,'bo-','LineWidth',2,'MarkerSize',msize)
load('diff_PUGSD_eps_0.01p4.mat')
semilogy(1:length(errresPUGSD),errresPUGSD,'r^-','LineWidth',2,'MarkerSize',msize)
load('diff_PUGSD_eps_0.001p4.mat')
semilogy(1:length(errresPUGSD),errresPUGSD,'k>-','LineWidth',2,'MarkerSize',msize)
legend('\epsilon_s = 5.10^{-1}','\epsilon_s = 10^{-1}','\epsilon_s = 10^{-2}','\epsilon_s = 10^{-3}')
set(gca,'FontSize',12)
xlabel('M')
ylabel('|| R_M ||')

%%
uM=uMPUGSD
% uM=spectral_decomposition(uM);
uref = uMPUGSD;
for i=1:getm(uM)
    errL2PUGSD(i) = norm(uref-truncate(uM,1:i))/norm(uref);
end
uM=uMPGSD
% uM=spectral_decomposition(uM);
uref = uMPGSD;
for i=1:getm(uM)
    errL2PGSD(i) = norm(uref-truncate(uM,1:i))/norm(uref);
end
%%



nbmodes=nbfun-2;
figure(44)
clf
semilogy(1:nbmodes,errL2PGSD(1:nbmodes),'b--s','LineWidth',2)
hold on
semilogy(1:nbmodes,errL2PUGSD(1:nbmodes),'r-o','LineWidth',2)
set(gca,'FontSize',fsize)
set(gca,'Xtick',[1:nbmodes])
xlim([.8,nbmodes])
xlabel('M')
ylabel('|| u_M - u_{ref} ||')
legend('algorithm 1','algorithm 2')
%%
if dim==1
    uq = quantile(upc,[0.1,0.9],3e4);
    figure(100)
    plotenveloppe(getcoord(getnode(S)),unfreevector(S,uq),'y')
    hold on
    figure(101)
    pdfplot(expand(uM(25)),'nbs',2e4,'b','ksdensity')
    hold on
    figure(102)
    pdfplot(expand(uM(40)),'nbs',2e4,'b','ksdensity')
    hold on
end

if dim==2
    rep1=getnodenextto(getnode(S),POINT([1,3/2]));
    rep1 = findddl(S,'u',rep1) ;
    rep2=getnodenextto(getnode(S),POINT([.5,1/3]));
    rep2 = findddl(S,'u',rep2) ;
    
    figure(101)
    pdfplot(expand(upc(rep1)),'nbs',2e4,'b','ksdensity')
    hold on
    figure(102)
    pdfplot(expand(upc(rep2)),'nbs',2e4,'b','ksdensity')
    hold on
end

%%
vi=1;
nc=3;
for up=[0,1]
    for ort = [0,1]
        for nbmodes = [12]
            
            nbmodes=min(nbmodes,getm(uM));
            nl=nbmodes/nc;
            
            
            fich = ['fig_diff'];
            if nc==3
                fich = [fich , '_3col'];
            end
            if up
                uM=uMPUGSD;
                fich = [fich , '_PUGSD'];
            else
                uM=uMPGSD;
                fich = [fich, '_PGSD'];
            end
            
            fich = [fich, '_' num2str(nbmodes) 'modes'];
            if ort
                stoU = double(getV(uM));[stoU,R] = qr(stoU,0);uM = PCRADIALMATRIX(stoU,size(uM),stol);
                fich = [fich , '_ortho'];
            end
            
            
            figure(16)
            clf
            
            if vi==2
                plotV(uM,1:nbmodes,S,'nl',nl,'nc',nc,'fact',.8,'manutext',{{0.6,-1},'U','FontSize',16},'EdgeColor','none','LineWidth',0)
            else
                plotV(uM,1:nbmodes,S,'nl',nl,'nc',nc,'surface','fact',1,'manutext',{{0.6,-1},'U','FontSize',16},'EdgeColor','none','LineWidth',0)
            end
            
            mysaveas('E:\REDACTION\matlab_results\Lshape_nonlin\', fich ,{'epsc2','jpeg'});
            
        end
    end
end
%%
vi=1
nc=2
nbmodes=8
uM=truncate(uMrefL2,1:nbmodes);
nbmodes=min(nbmodes,getm(uM));
nl=nbmodes/nc;

fich = ['fig_diff_PUGSD'];
fich = [fich, '_' num2str(nbmodes) 'modes'];

figure(16)
clf

if vi==2
    plotV(uM,1:nbmodes,S,'nl',nl,'nc',nc,'fact',.8,'manutext',{{0.6,-1},'U','FontSize',16},'EdgeColor','none','LineWidth',0)
else
    plotV(uM,1:nbmodes,S,'nl',nl,'nc',nc,'surface','fact',1,'manutext',{{0.6,-1},'U','FontSize',16},'EdgeColor','none','LineWidth',0)
end
mysaveas('C:\REDACTION\matlab_results\Lshape_nonlin\', fich ,{'epsc2','jpeg'});


%%
nbmodes=5;
fich = ['fig_diff'];

uMtest{1} = unfreevector(S,uMPGSD);
uMtest{2} = unfreevector(S,uMPUGSD);

myfunplot = @(x) abs(x);% log(abs(x));
figure(8)
clf
scanj=[1,2,4,7]
for j=scanj
    for k=1:2
        fullsubplot(length(scanj),2,2*(find(j==scanj)-1)+k)
        uM = uMtest{k};
        plot(myfunplot(((mean(uM)-mean(truncate(uM,1:j)))./max(abs(mean(uM))))),S);
        % axis on
        % if j==1;ax0=axis;cax0=caxis();else;axis(ax0);caxis(cax0);end
        % set(gca,'FontSize',12)
        
        text(1.4,.1,[ 'u' '_{' num2str(j) '}'],'FontSize',14);
        if k==1
            cax0=caxis;
            ax0=axis;
        else
            caxis(cax0)
            axis(ax0)
        end
        colormap('default')
        % colormap(flipud(gray))
        colorbar
        
    end
end

figure(9)
clf
for j=scanj
    for k=1:2
        fullsubplot(length(scanj),2,2*(find(j==scanj)-1)+k)
        uM = uMtest{k};
        plot(myfunplot((std(uM)-std(truncate(uM,1:j)))./max(abs(std(uM)))),S);
        % set(gca,'FontSize',12)
        % if j==1;ax0=axis;cax0=caxis();else;axis(ax0);caxis(cax0);end
        
        text(1.4,.1,[ 'u' '_{' num2str(j) '}'],'FontSize',14);
        if k==1
            cax0=caxis;
            ax0=axis;
        else
            caxis(cax0)
            axis(ax0)
        end
        colormap('default')
        colorbar
        
    end
    
end

%%
fich = ['fig_diff'];

uMtest{1} = uMPGSD;
uMtest{2} = uMPUGSD;

figure(9)
clf
figure(10)
clf
scanj=[2,6,10,16]
for j=scanj
    for k=1:2
        % subplot(length(scanj),2,2*(find(j==scanj)-1)+k)
        figure(8+k)
        subplot(length(scanj),1,find(j==scanj))
        uM = uMtest{k};
        RM = fun_res(truncate(uM,1:j));
        if k==1
            fprintf('PGSD')
            
        else
            fprintf('PUGSD')
            
        end
        norm(RM)
        RMmom2 = (std(RM).^2+mean(RM).^2);
        plot(RMmom2,S,'surface');
        % set(gca,'FontSize',12)
        % if j==1;ax0=axis;cax0=caxis();else;axis(ax0);caxis(cax0);end
        
        text(1.4,.1,[ 'u' '_{' num2str(j) '}'],'FontSize',14);
        if k==1
            cax0=caxis;
            ax0=axis;
        else
            caxis(cax0)
            axis(ax0)
        end
        colormap('default')
        colorbar
        
        pause(.5)
    end
    
end


%%

figure(10)
clf
for j=1:getm(uM)
    subplot(5,2,j)
    ax = full(mean(getL(uM,j))+3*[-1,1]*std(getL(uM,j)));
    pdfplot(getL(uM,j),'nbs',1e5,'ksdensity','npts',200,'axis',ax);
    % if j==1;ax0=axis;cax0=caxis();else;axis(ax0);caxis(cax0);end
    colormap(flipud(gray))
    pause(.01)
end


%%


%% ARNOLDI
nbfun = 12;
orthocrit=1e-11;
uM = PCRADIALMATRIX([size(u)],PC);
uM=full(uM);
testerror=0;
nbreac = 1;
nbiterinterne=3;
clear errresAGSD
for i=0:nbreac
    fprintf(' ARNOLDI reac #%d \n',i)
    lam = allones(PC);lam=lam/norm(lam);
    W = double(getV(uM));
    for k=1:nbfun-getm(uM)
        %    for j=1:nbiterinterne
        U = fun_FM(lam,uM);
        
        U=U/norm(U);
        U = U-W*(W'*U);
        fprintf(' ARNOLDI iter #%d : residu orthogonalisation = %.3d\n',k,norm(U))
        H = norm(U);
        U = U/H;
        lam = fun_GM(U,uM);
        
        lam = lam/norm(lam);
        %    end
        if H<orthocrit
            fprintf(' break\n')
            break
        end
        
        
        W = [W,U];
        stol = [getL(uM);lam];
        [W,R] = qr(W,0);
        stol = R*stol;
        
        if testerror && k<nbfun-getm(uM)
            uMtemp = PCRADIALMATRIX(W,size(u),fun_G(W,stol));
            errresAGSD(getm(uMtemp)) =  fun_reserror(uMtemp);
            fprintf(' GSD ARNOLDI : number of functions = %d , erreur residu %.3d\n',getm(uMtemp),errresAGSD(getm(uMtemp)));
        end
    end
    
    
    if testerror && k~=nbfun-getm(uM)
        uM = uMtemp;
    else
        uM = PCRADIALMATRIX(W,size(u),fun_G(W,[]));
        errresAGSD(getm(uM)) = fun_reserror(uM);
        fprintf(' GSD ARNOLDI : number of functions = %d , erreur residu %.3d\n',getm(uM),errresAGSD(getm(uM)));
        
    end
    
    if getm(uM)>=nbfun
        break
    end
    
end



%% affichage des quantites

qplotquan = unfreevector(S,uM);
% qplotquan = qradreuse-expect(qradreuse);
% qplotquan = orthogonalizeL(qplotquan);
% qplotquan=normsort(qplotquan);
choixdroite=2;
switch choixdroite
    case 1
        D10 = LINE(POINT([1,0]),POINT([1,2]));
        [rep,P10]=ispointin(D10,POINT(getnode(S)));
        t=double(getcoord(P10));t=t(1,2,:);t=t(:);t=sort(t);
        P10 = POINT([ones(length(t),1),t]);
    case 2
        D10 = LINE(POINT([0,0]),POINT([0,2]));
        [rep,P10]=ispointin(D10,POINT(getnode(S)));
        t=double(getcoord(P10));t=t(1,2,:);t=t(:);t=sort(t);
        P10 = POINT([zeros(length(t),1),t]);
end
qdroitequan=eval_sol(S,qplotquan,P10,'u');


m = getm(qplotquan);
nl = ceil(sqrt(m));nc = ceil(m/nl);

Ymean = expect(qdroitequan);
fillcolor={'k','b','r','m','g','c','y'};
fillcolor=fillcolor(end:-1:1);
nr = 10000;
quan = [0.01,0.99];
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



figure(206*choixdroite)
clf
hold on
leg = {};
V = getmode(qdroitequan,1:nbmodes);
Yref = quantile(V,quan,nr);
scanmode = [1,2,3,4];
nl = ceil(sqrt(length(scanmode)));nc = ceil(length(scanmode)/nl);
for k=1:length(scanmode)
    subplot(nl,nc,k)
    i=scanmode(k);
    V = getmode(qdroitequan,1:i);
    Y = quantile(V,quan,nr);
    plotenveloppe(t,Yref,'w')
    hold on
    plotenveloppe(t,Y,'y');
    leg = [{'reference'} , {['mode 1 to ' num2str(i)]}];
    legend(leg{:})
end


%%
uM = unfreevector(S,uM)

rep1=getnodenextto(getnode(S),POINT([1.5,1.5]));
rep1 = findddl(S,'u',rep1) ;
rep2=getnodenextto(getnode(S),POINT([.5,0.1]));
rep2 = findddl(S,'u',rep2) ;
options={'nbs',2e5,'npts',30,'ksdensity','LineWidth',2}
rep=rep1
figure(102)
clf
pdfplot(expand(uM(rep)),'k-',options{:})
ax=axis;
pause(.1)
hold on
pdfplot(expand(getmode(uM(rep),1:1)),'r--',options{:})
axis(ax);
pause(.1)
pdfplot(expand(getmode(uM(rep),1:2)),'b:',options{:})
axis(ax);
pause(.1)
pdfplot(expand(getmode(uM(rep),1:3)),'g-.',options{:})
axis(ax);
pause(.1)
pdfplot(expand(getmode(uM(rep),1:5)),'m--',options{:})
axis(ax);
set(gca,'FontSize',14)
legend('reference','M = 1','M = 2','M = 3','M = 5')
rep=rep2

figure(103)
clf
pdfplot(expand(uM(rep)),'k-',options{:})
ax=axis;
pause(.1)
hold on
pdfplot(expand(getmode(uM(rep),1:1)),'r--',options{:})
axis(ax);
pause(.1)
pdfplot(expand(getmode(uM(rep),1:2)),'b:',options{:})
axis(ax);
pause(.1)
pdfplot(expand(getmode(uM(rep),1:3)),'g-.',options{:})
axis(ax);
pause(.1)
pdfplot(expand(getmode(uM(rep),1:5)),'m--',options{:})
axis(ax);
set(gca,'FontSize',14)
legend('reference','M = 1','M = 2','M = 3','M = 5')




