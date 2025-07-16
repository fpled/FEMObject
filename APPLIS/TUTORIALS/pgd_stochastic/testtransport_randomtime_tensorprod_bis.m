
vx = @(x) (x(:,2)-1/2);
vy = @(x) (1/2-x(:,1));

S = gmsh2femobject(2,[getfemobjectoptions('path') 'APPLIS/TUTORIALS/pgd_stochastic/example_transport2.geo'],2);
xnode = getcoord(getnode(S));
VX = vx(xnode);
VY = vy(xnode);

mat = FOUR_ISOT('k',1,'b',...
    {{FENODEFIELD(VX),FENODEFIELD(VY)}},'c',1);
S = setmaterial(S,mat);
S=final(S);
S=addcl(S,[],'T');
Mx = calc_freematrix(S,@mass);
Adiff = calc_freematrix(S,@diff);
Aadv = 250*calc_freematrix(S,@adv);
Areac = 10*calc_freematrix(S,@mass);
fx = bodyload(keepgroupelem(S,2),[],'QN',100);
foutput = bodyload(keepgroupelem(S,3),[],'QN',1);

T = TIMEMODEL(0,.03,80);
N = DGTIMESOLVER(T,0);
Dt = getDmatrix(N,'basic');
Mt = getMmatrix(N);
Dtdiff = getDmatrix(N);

ftone = full(double(one(N))');
ft = getMmatrix(N)*ftone;

figure(3)
clf
plot(S)

%%

utdet = dsolve(N,fx*one(N),Mx,Adiff+Aadv+Areac);

figure(8)
clf
mova = evol(utdet,S,'surface','rescale','z')

figure(9)
clf
plot(foutput'*utdet);
integrate(foutput'*utdet)
%%
tp = 0;

casspeed = 1;
%
m1=1;
m2=m1;
m3=m1;
m4=0;
%m5=2;
p=4;
sigma1=.2;sigma2=.2;sigma3=.2;%sigma4=.2;%sigma5=.3;
RV1 = setnumber(RVUNIFORM(1-sigma1,1+sigma1),1);
RV2 = setnumber(RVUNIFORM(1-sigma2,1+sigma2),2);
if casspeed==1
    RV3 = setnumber(RVUNIFORM(1-sigma3,1+sigma3),3);
    options={};
elseif casspeed==2
    RV3 = setnumber(RVUNIFORM(-3,1),3);
    options = {'fedim',3,'femesh',{4}};
end
RV = RANDVARS(RV1,RV2,RV3);

block = 0;

if block==0
    typebase = 1;
    [X,PC] = PCTPMODEL(RV,'order',[p,p,p],'pcg','typebase',typebase,options{:});
    a1 = X{1} ; a2=X{2} ; a3=X{3} ; a4=one(PC);
elseif block==1
    typebase = 2;
    groups = {[1,2],3};
    [X,PC] = PCTPMODEL(RV,'order',[p,p,p],'pcg','groups',groups,'typebase',typebase,options{:});
    a1 = X{1} ; a2=X{2} ; a3=X{3} ; a4=one(PC);

elseif block==2
    typebase = 2;
    groups = {1:3};
    [X,PC] = PCTPMODEL(RV,'order',[p,p,p],'pcg','groups',groups,'typebase',typebase,'fedim',3,'femesh',{4});
    a1 = X{1} ; a2=X{2} ; a3=X{3} ; a4=one(PC);

end
%
M = full(one(PC))*Mx;
A = a1*Adiff + a3*Aadv + a2*Areac;
A = calc_ximasse(A);
b = a4*fx;

if block==0
    Mone = get_ximasse(calc_ximasse(one(PC)));
    Msep = SEPMATRIX([{Mx},Mone]);
    Asep = SEPMATRIX(getnbgroups(PC)+1);
%Asep = Asep + SEPMATRIX([{Adiff+Aadv+Areac},Mone]);
    bsep = SEPMATRIX([{fx*ftone'},getphi(one(PC))]);

    Mi = get_ximasse(X{1});
    Asep = Asep + SEPMATRIX([{Adiff},Mi]);
    Mi = get_ximasse(X{3});
    Asep = Asep + SEPMATRIX([{Aadv},Mi]);
    Mi = get_ximasse(X{2});
    Asep = Asep + SEPMATRIX([{Areac},Mi]);

%for i=1:m4
%bsep = bsep + SEPMATRIX([{fx*ftone'},getphi(X{m1+m2+m3+i})]);
%end

else
    Mone = get_ximasse(calc_ximasse(one(PC)));
    Msep = SEPMATRIX([{Mx},Mone]);
    a1M = get_ximasse(calc_ximasse(a1));
    a2M = get_ximasse(calc_ximasse(a2));
    a3M = get_ximasse(calc_ximasse(a3));
    Asep = SEPMATRIX([{Adiff},a1M]) + ...
        SEPMATRIX([{Aadv},a3M]) + ...
        SEPMATRIX([{Areac},a2M]);
    bsep = SEPMATRIX([{fx*ftone'},getphi(a4)]);

end

%%
ur = dsolve(N,random(b)*one(N),random(M),random(A));
figure(4)
clf
evol(ur,S,'surface','rescale','z');


%% monte carlo

nbs = 1e4;
omc = [];
omctime=zeros(0,length(N));
%%
fich = 'testtransport_random_bis_m3_sigma_montecarlo';
tic
for kk=1:nbs
    pourcentage(kk,nbs,100)
    utemp = dsolve(N,random(b)*one(N),random(M),random(A));
    utemp = foutput'*double(utemp);
    omc = [omc,utemp*ft];
    omctime=[omctime;utemp];
    if mod(length(omc),500)==0
        save(fich,'omc','omctime')
    end
end
toc
save(fich,'omc','omctime')

%
%%
figure(2)
hold on
pdfsampleplot(omc,'k','LineWidth',2)%,'npts',20)

%%
scanpfval = [4:.2:5.5]*10^-4;
pfmc = zeros(1,length(scanpfval));
for i=1:length(scanpfval)
    pfmc(i) = length(find(omc>scanpfval(i)))/length(omc);
end



%%
solref=[];
maxorder = 20;
maxiter=5;
update = 1;
alphaupdate = 0;
updatedim  =  2 :getdim(Asep);
adjoint =0;
tol=1e-9;
updatestep=1;
itercritupdate=1e-4;
col = 'm--^';%col = 'b-s';%col='b--d'
%col = 'k-o'
errorindicator = 'none';
%errorindicator = 'stagnation';

solver = SEPSOLVER(getdim(Asep),'tol',tol,'alphaupdate',alphaupdate,...
    'update',update,'updatedim',updatedim,...
    'maxorder',maxorder,'maxiter',maxiter,'reference',solref,...
    'adjoint',adjoint,'updatestep',updatestep,...
    'itercritupdate',itercritupdate,...
    'errorindicator',errorindicator,'itercrit',1e-3);
tic
[usep,resultsep] = dsolve(Msep,Asep,bsep,N,solver);
toc
figure(17);semilogy(resultsep.error,col,'LineWidth',2);
hold on
x0 = xlim;
set(gca,'FontSize',14)
xlim([1,max(x0(2),length(resultsep.error))]);
useppc = PCTPMATRIX(usep,PC);
oseptime = foutput'*useppc;
osep = foutput'*useppc*ft;

eval(['useppc' num2str(maxorder) '=useppc;']);
%%

figure(2)
hold on
pdfplot(osep,col,'nbs',1e5,'LineWidth',2,'npts',20)
set(gca,'FontSize',14)

%%
figure(4)
clf
[urr,xi]=random(useppc);
urr = TIMEMATRIX(urr,N);
randomeval(a3,xi)
evol(urr,S,'surface','rescale','z');

%%
ur = dsolve(N,randomeval(b,xi)*one(N),randomeval(M,xi),randomeval(A,xi));
figure(4)
clf
evol(ur,S,'surface','rescale','z');



%%

xi = [zeros(100,1),linspace(-1,1,100)',zeros(100,1)];
osepxi = randomeval(osep,xi)
figure(4)
plot(xi(:,2),osepxi)


%% COMPARAISON SAMPLES AVEC REFERENCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
N = setevolparam(N,'plotstep',10);

[ur,xi] = random(useppc);
ur = TIMEMATRIX(ur,N);
Ar = randomeval(A,xi);
Mr = randomeval(M,xi);
br = randomeval(b,xi);
urref = dsolve(N,br*one(N),Mr,Ar);

figure(3)
%evol(TIMEMATRIX(ur,N),S,'surface','rescale','z')
evolV(N,{urref,ur},1:2,S,'surface','rescale','z')

%%
useptest = {useppc1,useppc2,useppc4,useppc8,useppc12};
%useptest={useppc};
normerr=cell(1,length(useptest));
normref=cell(1,length(useptest));
nbs=50;
for kk=1:nbs
    pourcentage(kk,nbs,nbs)
    [ur,xi] = random(useptest{1});
    ur = TIMEMATRIX(ur,N);
    Ar = randomeval(A,xi);
    Mr = randomeval(M,xi);
    br = randomeval(b,xi);
    urref = dsolve(N,br*one(N),Mr,Ar);
    normerr{1}(end+1) = norm(double(ur)-double(urref));
    normref{1}(end+1) = norm(double(urref));
    for ii=2:length(useptest)
        ur = randomeval(useptest{ii},xi);
        normerr{ii}(end+1) = norm(double(ur)-double(urref));
        normref{ii}(end+1) = norm(double(urref));
    end
end

errsup = zeros(1,length(normref));
errmin = zeros(1,length(normref));
errL2 = zeros(1,length(normref));
for ii=1:length(useptest)
    errsup(ii) = max(normerr{ii}./normref{ii});
    errmin(ii) = min(normerr{ii}./normref{ii});
    errL2(ii) = sqrt(sum(normerr{ii}.^2)/sum(normref{ii}.^2));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PDF DE L'OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2)
clf
pdfsampleplot(omc,'k','LineWidth',2)%,'npts',20)
hold on
oseptest = {useppc1,useppc2,useppc4,useppc8,useppc12};
for j=1:length(oseptest)
    oseptest{j} =  foutput'*oseptest{j}*ft  ;
end
scancol = {'m-^','b-d','r-o','c-v'};

for j=1:length(scancol)
    pdfplot(oseptest{j},scancol{j},'nbs',1e5,'LineWidth',2,'npts',20)
    hold on
    set(gca,'FontSize',14)
    pause(.01)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AFFICHAGE DES MODES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = setevolparam(N,'plotstep',2);
N = setevolparam(N,'plottime',0);
figure(4)
clf
usepplot = usep.F(:,1);
usepplot{2}=-usepplot{2};
mova = evolV(N,usepplot,1:length(usepplot),S,'surface');
fich = './illustreseparation/testtransport_time_1block_';
fich = [fich 'modes1_to_' num2str(length(usepplot)) ];
movie2avi(mova,fich)

%%
usepplot = usep.F(:,1);
usepplot{2}=-usepplot{2};
N = setevolparam(N,'plotstep',2);
N = setevolparam(N,'plottime',0);
figure(5)
clf
for ii=1:8
    mova = evol(N,TIMEMATRIX(usepplot{ii},N),S,'surface','rescale','z');
%mova = evol(N,TIMEMATRIX(usepplot{ii},N),S);
    fich = './illustreseparation/testtransport_time_1block_';
    fich = [fich 'surf_mode' num2str(ii) ];
    movie2avi(mova,fich)
end

%%
repttplot = 10:20:80
for ii=1:8
    utemp = usepplot{ii}(:,repttplot);
    cax = [min(min(utemp)),max(max(utemp))];
    figure(6)
    clf
    for k=1:length(repttplot)
        subplot(1,length(repttplot),k)
        plot(utemp(:,k),S)
        caxis(cax);
    end
    fich = './illustreseparation/testtransport_time_1block_';
    fich = [fich 'mode' num2str(ii) '_cliches' ];
    myprint('',fich,'epsc2')
end


%% MIX GSD

solref=[];
col='g-s'

m=10;
update = 1;
updatel=1;
ml=1;
updatelu=1;
mlu = 30;
itercrit = 1e-2;
dimgsd = 1;
updatedim = 2:getdim(Asep);
errorindicator = 'none';
sv = SEPSOLVER(getdim(Asep),'tol',1e-3,...
    'update',update,'maxorder',m,'maxiter',5,'reference',solref,...
    'itercrit',itercrit,'updatestep',1,'errorindicator',errorindicator,...
    'type','arnoldi','restart',2,'fullupdate',true,...
    'orthocrit',1e-10);

svl = SEPSOLVER(getdim(Asep),'tol',1e-2,...
    'update',updatel,'maxorder',ml,'maxiter',5,...
    'updatedim',updatedim,'display',false,'inittype','rand',...
    'itercrit',itercrit);
svlu = svl;
svlu = setparam(svlu,'maxorder',mlu,'maxiter',5,...
    'update',updatelu,'updatedim',updatedim,'display',true,...
    'itercrit',itercrit,'tol',1e-4,'updatestep',1,'errorindicator','none');

[usep,resultgsd,W] = dsolve_gsd_arnoldi(Msep,Asep,bsep,N,sv,svl,svlu);

figure(17);semilogy(resultgsd.error,col);
hold on

useppc = PCTPMATRIX(usep,PC);
eval(['useppc' num2str(m) '=useppc;']);
oseptime = foutput'*useppc;
osep = foutput'*useppc*ft;


