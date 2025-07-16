
multi_inclusions
S = createddlnode(S,DDL('u'));
S = addcl(S,[],'u');

%%
fx = LINFORM(0,1);
fx = setselgroup(fx,getnbgroupelem(S));
fx = fx{S}(:);
out = fx;

Ax = cell(1,getnbgroupelem(S));
for i=1:getnbgroupelem(S)
    Ax{i} = BILINFORM(1,1);
    Ax{i} = setselgroup(Ax{i},i);
    Ax{i} = Ax{i}{S}(:,:);
end

RV = RANDVARS(RVLOGUNIFORM(1e-1,2),getnbgroupelem(S)-2);
RV{1} = RVLOGUNIFORM(1e-2,10);
p=10;
[X,PC] = PCTPMODEL(RV,'order',p,'pcg','typebase',2);

A = (Ax{1}+Ax{end})*one(PC);
for i=1:getnbgroupelem(S)-2
    A = A + Ax{i+1}*X{i};
end
b = fx*one(PC);

Asep = SEPMATRIX(A);
bsep = SEPVECTOR(b);

%%  calcul progressif

PGD = SEPSOLVER(getdim(Asep),'maxorder',30,'tol',1e-5,'display',true,...
    'maxiter',10,'updatedim',2:getdim(Asep),'update',1);
[usep,resultref] = solve(Asep,bsep,PGD);
useppc = PCTPMATRIX(usep,PC,2:getdim(usep));
oseppc=out'*useppc;
osep = mtimes(out',usep,1);
oref=osep;
uref = usep;
urefpc = useppc;
orefpc = out'*useppc;
sobref = sobol_indices(orefpc);
figure(12)
clf
semilogy(resultref.error)

save('save_test_separation_multinclusions','uref','urefpc','orefpc','sobref','oref')


%% convergence des indices de sobol
rrr=[];
for kk=1:getm(oseppc)
    rrr = [rrr;sobol_indices(truncate(oseppc,1:kk))];
end
figure(45)
plot(rrr)
%% calcul sur bases r�duites
WAsep = mtimes(W',Asep,1);
WAsepW=mtimes(WAsep,W,1);
Wbsep=mtimes(W',bsep,1);

PGD = SEPSOLVER(getdim(Asep),'maxorder',50,'tol',1e-5,'display',true,'maxiter',10);
[usep,result] = solve(WAsepW,Wbsep,PGD);
usep = mtimes(W,usep,1);
useppc = PCTPMATRIX(usep,PC,2:getdim(usep));
oseppc=out'*useppc;

figure(12)
clf
semilogy(result.error)

%% resolution avec 2 niveaux de tensor product

PGD = SEPSOLVER(getdim(Asep),'maxorder',20,'tol',1e-5,'display',true,'maxiter',10,'type','alterne',...
    'update',1,'updatestep',3);
PGDl =  SEPSOLVER(getdim(Asep),'maxorder',5,'tol',1e-2,'display',false,'maxiter',10);
PGDL =  SEPSOLVER(getdim(Asep),'maxorder',50,'tol',1e-3,'display',true,'maxiter',10);
[usep,result] = solve_gsd(1,Asep,bsep,PGD,PGDl,PGDL);
useppc = PCTPMATRIX(usep,PC,2:getdim(usep));
oseppc=out'*useppc;

%% simulations de montecarlo

stoos=[];
stors=[];
stooseps=[];

%% sur modele full
nbsamples = 2000
tic
for kk=1:nbsamples
    pourcentage(kk,nbsamples,100);
    %grad=1;
    %bi = BILINFORM(1,1);
    [As,rs]=random(A);
    stors=[stors;rs];
    bs = randomeval(b,rs);
    us = As\bs;
    os = out'*us;
    stoos = full([stoos,full(os)]);
end
toc
Q = length(stoos);
EX = sum(stoos)/Q;
VX = sum(stoos.*stoos)/Q - EX^2;
%% sur modele reduit
sr = 5;
V = zeros(size(A,1),sr);
for kk=1:sr
    [As,rs]=random(A);
    stors=[stors;rs];
    bs = randomeval(b,rs);
    V(:,kk) = As\bs;
end
V = orth(V);
Ar = V'*A*V;br = V'*b;
outr = V'*out;
nbsamples = 100;
tic
for kk=1:nbsamples
    pourcentage(kk,nbsamples,10);
    %grad=1;
    %bi = BILINFORM(1,1);
    [As,rs]=random(Ar);
    stors=[stors;rs];
    bs = randomeval(br,rs);
    us = As\bs;
    os = outr'*us;
    stoos = full([stoos,full(os)]);
end
toc
    

%%
stosob = [];
for var=1:8
    stoosvar = [];
    storsvar = [];

    for kk=1:length(stoos)
        pourcentage(kk,length(stoos))
        bi = BILINFORM(1,1);
        [As,rs]=random(A);
        rs(var)=stors(kk,var);
        storsvar=[storsvar;rs];
        As =randomeval(A,rs);
        bs = randomeval(b,rs);
        us = As\bs;
        os = out'*us;
        stoosvar = full([stoosvar,os]);
    end
%
    Q = length(stoos);
    EX = sum(stoos)/Q;
    VX = sum(stoos.*stoos)/Q - EX^2;
    Vi = sum(stoos.*stoosvar)/Q - EX^2;
    stosob(var) = Vi/VX
end

%%
bi = BILINFORM(1,1);
[As,rs]=random(A);bs = randomeval(b,rs);
us = As\bs;
os = out'*us;
grad=1;
figure(1)
clf
if grad>=1
    plot_sol(S,us,'epsilon',grad);
else
    plot(us,S);
end
cax = caxis;
figure(2)
clf
useps=randomeval(useppc,rs);
oseps = out'*useps;
if grad>=1
    plot_sol(S,us,'epsilon',grad);
else
    plot(useps,S)
end
caxis(cax)
norm(us-useps)/norm(us)
sqrt(bi{S}(us-useps,us-useps)/bi{S}(us,us))
full(abs(oseps-os)/abs(os))


figure(56)
clf
plot(std(useppc),S)

%%
figure(45)
clf
plot(stoos,stooseps,'*')
hold on
plot(xlim,xlim,'k-')


%%
Ws =[];
for k=1:100
    [As,rs]=random(A);
    bs = randomeval(b,rs);
    us = As\bs;
    Ws = [Ws,us];
end
[W,Dtemp,Vtemp,err]=svdtruncate(full(Ws),1e-3);


%% affichage des indices de sobol

figure(3)
clf
plot(S)
for i=1:length(sobref)
    plot(S,'selgroup',i+1,'Color',sobref(i))
end
colormap(cool)
colorbar
