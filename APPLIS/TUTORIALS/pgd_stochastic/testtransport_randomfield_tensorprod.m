%%%%%%%%%%%%%%%%%%%%%%%%
% EXEMPLE 2
% ARTICLE NOUY., HIGH DIMENSIONAL STOCHASTIC, ARCHIVES COMP, 2010
%%%%%%%%%%%%%%%%%%%%%%%%
%
% time dimension is considered as a general dimension
%%%%%%%%%%%%%%%%%%%%%%%%%


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
ft = getMmatrix(N)*full(double(one(N))');

figure(3)
clf
plot(S)


%%
[Skl,Xkl,Ykl] = mesh(DOMAIN(2),40,40);

mkl=100;
[Vkl,Dkl] = KL(RFGAUSSIANEXP2(1,0.2,.01),mkl,Skl,'modes');
Vkl=Vkl*Dkl;
Dkl = diag(Dkl);
xnode = getcoord(getnode(S));
V=zeros(S.nbnode,mkl+1);
for kk=1:mkl+1
    V(:,kk) = interp2(Xkl,Ykl,reshape(Vkl(:,kk),size(Xkl)),xnode(:,1),xnode(:,2));
end
figure(56)
clf
multiplot(mat2cell(V,size(V,1),ones(1,mkl+1)),1:30,S)
%eval_sol(V)

figure(45)
clf
semilogy((Dkl(2:end)),'*','MarkerSize',8)
set(gca,'FontSize',16)
ylabel('(\sigma_i)^{1/2}')
xlabel('i')
fich = './illustreseparation/testtransport_randomfield_m40_diffusion_randfield_spectrum';
myprint('',fich,'eps')
%%
onlykl=1;
RV = RANDVARS(RVUNIFORM(-1,1),mkl);
sg = 1;
p = 4;
groups={};
nbgroups=0;
alldim = 1:mkl;
while ~isempty(alldim)
    groups{nbgroups+1}=alldim(1:min(sg,end));
    alldim = setdiff(alldim,alldim(1:min(sg,end)));
    nbgroups=nbgroups+1;
end

if ~onlykl
    groups = [groups,{mkl+1:mkl+2}];
end

if mkl==1
    groups={[groups{:}]};
end
[X,PC] = PCTPMODEL(RV,'order',p,'pcg','groups',groups);

kmat = V(:,1)*one(PC);
for i=1:mkl
    kmat = kmat+V(:,i+1)*X{i};
end
if onlykl
    cmat=one(PC);
    rmat=one(PC);
else
    cmat = X{mkl+1};
    rmat = X{mkl+2};
end

cmat=calc_ximasse(cmat);
rmat=calc_ximasse(rmat);


%%

A = Aadv*cmat + Areac*rmat ;
A = A+Adiff*one(PC);
for i=1:mkl
    pourcentage(i,mkl)
    a = BILINFORM(1,1,V(:,i+1),0);
    A = A+ a{S}(:,:)*X{i};
end
M = full(one(PC))*Mx;
M = calc_ximasse(M);
A = calc_ximasse(A);
b = one(PC)*fx;

%%

Mone = get_ximasse(calc_ximasse(one(PC)));
Msep = SEPMATRIX([{Mx},Mone]);
Asep = SEPMATRIX([{Adiff},Mone]);
for i=1:mkl
    pourcentage(i,length(X)-2)
    Mi = get_ximasse(X{i});
    a = BILINFORM(1,1,V(:,i+1),0);
    Adiffi = a{S}(:,:);
    Asep = Asep + SEPMATRIX([{Adiffi},Mi]);
end
Mc = get_ximasse(cmat);
Asep = Asep + SEPMATRIX([{Aadv},Mc]);
Mr = get_ximasse(rmat);
Asep = Asep + SEPMATRIX([{Areac},Mr]);
bsep = SEPMATRIX([{fx},getphi(one(PC))]);



%%
adv1 = 1;
areac1 = 1;
fx = bodyload(keepgroupelem(S,2),[],'QN',100);
umoy = (Adiff+adv1*Aadv+areac1*Areac)\fx;

kmatr = random(kmat);
cmatr = adv1*random(cmat);
rmatr = areac1*random(rmat);

figure(66)
clf
plot(kmatr,S,'surface')
colorbar

a = BILINFORM(1,1,kmatr,0);
Adiffr = a{S}(:,:);
ur = (Adiffr+cmatr*Aadv+rmatr*Areac)\fx;
figure(67)
clf
subplot(1,2,1)
plot(ur,S)
cax=caxis;
colorbar
subplot(1,2,2)
plot((ur-umoy),S)
%caxis(cax)
colorbar

figure(68)
clf
plot(ur,S)


%%
ursto = zeros(size(ur,1),0);
%%
nbs=300;
for kkk=1:nbs
    pourcentage(kkk,nbs,min(nbs,100))
    kmatr = random(kmat);
    cmatr = random(cmat);
    rmatr = random(rmat);
    if min(kmatr)<0.1
        warning('trop faible')
    else
        a = BILINFORM(1,1,kmatr,0);
        Adiffr = a{S}(:,:);
        ur = (Adiffr+cmatr*Aadv+rmatr*Areac)\fx;

        ursto = [ursto,ur];
    end
end


%% monte carlo

istime=0;

nbs = 1e5;
omc = [];
omctime=zeros(0,length(N));

%%

for kk=1:nbs
    pourcentage(kk,nbs,100)
    utemp = random(A)\random(b);
    omc = [omc,full(foutput'*utemp)];
end

save('test_transport_randfield_omc_statio','omc')


%%
figure(2)
hold on
pdfsampleplot(omc,'k','LineWidth',2,'npts',20)


%%
solref=[];
maxorder = 15;
maxiter=5;
update = 1;
alphaupdate = 0;
updatedim = 2:getdim(Asep);
adjoint =0;
tol=1e-9;
updatestep=1;
itercritupdate=1e-1;
col = 'r-*';
errorindicator = 'none';
%errorindicator = 'residual';

solver = SEPSOLVER(getdim(Asep),'tol',tol,'alphaupdate',alphaupdate,...
    'update',update,'updatedim',updatedim,...
    'maxorder',maxorder,'maxiter',maxiter,'reference',solref,...
    'adjoint',adjoint,'updatestep',updatestep,...
    'itercritupdate',itercritupdate,...
    'errorindicator',errorindicator,'itercrit',1e-3);
tic
[usep,resultsep] = solve(Asep,bsep,solver);
toc
figure(17);semilogy(resultsep.error,col);
hold on

useppc = PCTPMATRIX(usep,PC);

if adjoint && update
    upseppcadjointupdate=useppc;
elseif adjoint
    upseppcadjoint=useppc;
elseif update
    upseppcupdate=useppc;
end


osep = foutput'*useppc;

useptest = useppc;


%%
figure(2)
hold on
pdfplot(osep,col,'nbs',5e4,'LineWidth',2)


%%
errtest = [];

%% CALCUL ERREUR POUR STATIO
nsol = [];
nerr = [];
%%
nbserr=20;
for kk=1:nbserr
    pourcentage(kk,nbserr)
    [ur,xi] = random(useptest);
    Ar = randomeval(A,xi);
    br = randomeval(b,xi);
    urref = Ar\br;
    nsol=[nsol,norm(double(urref))^2];
    nerr=[nerr,norm(ur-double(urref))^2];
    errtest = sqrt(nerr./nsol);
    errL2 = sqrt(mean(nerr)/mean(nsol))
end
errL2

figure(4)
clf
subplot(1,3,1)
plot(urref,S,'surface');
caxis([0,0.15])
colorbar
view(2)
subplot(1,3,2)
plot(ur,S,'surface');
caxis([0,0.15])
colorbar
view(2)
subplot(1,3,3)
plot((ur-urref)./max(abs(urref)),S,'surface');
view(2)
colorbar


%% AFFICHAGE

[ur,xi] = random(useptest);
Ar = randomeval(A,xi);
br = randomeval(b,xi);
urref = Ar\br;
kmatr=randomeval(kmat,xi);

figure(4);clf;subplot(1,3,1);
plot(urref-umoy,S,'surface');cax=caxis;colorbar
view(2)
subplot(1,3,2)
plot(ur-umoy,S,'surface');
caxis(cax);colorbar;view(2)
subplot(1,3,3)
plot((ur-urref)./max(abs(urref)),S,'surface');;view(2);colorbar




%%
solref=[];
col='g-s'

m=20;
update = 1;
updatel=1;
ml=3;
updatelu=1;
mlu = 25;
itercrit = 1e-2;
dimgsd = 1;
updatedim = 1:getdim(Asep);
errorindicator = 'none';
sv = SEPSOLVER(getdim(Asep),'tol',1e-6,...
    'update',update,'maxorder',m,'maxiter',5,'reference',solref,...
    'itercrit',itercrit,'updatestep',1,'errorindicator',errorindicator,...
    'type','arnoldi','restart',2,'fullupdate',true);

svl = SEPSOLVER(getdim(Asep),'tol',1e-2,...
    'update',updatel,'maxorder',ml,'maxiter',5,...
    'updatedim',updatedim,'display',false,'inittype','rand',...
    'itercrit',itercrit);
svlu = svl;
svlu = setparam(svlu,'maxorder',mlu,'maxiter',5,...
    'update',updatelu,'updatedim',updatedim,'display',true,...
    'itercrit',itercrit,'tol',5e-4,'updatestep',2,'errorindicator','none',...
    'adjoint',0);

tic
[usep,resultgsd,W,Lambda] = solve_gsd(dimgsd,Asep,bsep,sv,svl,svlu);
toc

useppc = PCTPMATRIX(usep,PC);
Lambdapc = PCTPMATRIX(Lambda,PC);
osep = foutput'*useppc;
useptest=useppc;

figure(17);semilogy(resultgsd.error,col);
hold on

eval(['useppc' num2str(m) '=useppc;']);
eval(['W' num2str(m) '=W;']);

save('test_transport_randfield_solution','useppc','W')


%% affichage des modes E(\mu\lambda^2) succesifs dans la resolution arnoldi

stolambda = resultgsd.stolambda;
stomu = cell(1,length(stolambda));
for ii=1:length(stolambda)
    stolambda{ii}=PCTPMATRIX(stolambda{ii},PC);
%stolambda{ii}=stolambda{ii}*(1/norm(stolambda{ii}));
    stomu{ii} = expect(kmat,stolambda{ii},stolambda{ii});
    stomu{ii} = stomu{ii}/expectmtimes(stolambda{ii},stolambda{ii});
end
%%
for ii=1:1
    figure(ii)
    clf
    plot(stomu{ii},S,'surface')
    view(2)
    colorbar
end
%%
figure(2)
hold on
pdfplot(ogsd,'b-','nbs',5e4,'LineWidth',2,'ksdensity')

%%

scanval = [3.5:.2:5]*10^-4;
pfsep = zeros(1,length(scanval));
pfmc = zeros(1,length(scanval));
pfgsd = zeros(1,length(scanval));
for i=1:length(scanval)
    pfmc(i) = length(find(omc>scanval(i)))/length(omc);
    pfsep(i) = length(find(osepr>scanval(i)))/length(osepr);
    pfgsd(i) = length(find(ogsdr>scanval(i)))/length(ogsdr);
end

figure(4)
clf
semilogy(scanval,pfmc,'k')
hold on
semilogy(scanval,pfsep,'b')
hold on
semilogy(scanval,pfgsd,'g')


