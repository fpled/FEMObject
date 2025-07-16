%% Definition du probleme
r=10;
P = POINT([ 0,0 ; 1,0 ; 0,2 ; 1,2 ; 1,1 ; 2,1 ; 2,2 ]);
S1 = mesh(DOMAIN(2,P(1),P(4)),r,r*2);
S2 = mesh(DOMAIN(2,P(5),P(7)),r,r);
S = union(S1,S2);
S = concatgroupelem(S);
RV=RANDVARS();

RV{1}=RVUNIFORM(0.8,1.2);
RV{2}=RVUNIFORM(2.4,3);
RV{3}=RVUNIFORM(2.4,3);
RV{4}=RVNORMAL(1,0.3);

X = PCMODEL(RV,'order',5,RANDPOLYS(RV));


mat=FOUR_ISOT('k',X{1},'b',[X{2};X{3}]);
S = setmaterial(S,mat);
S=final(S);
S=addcl(S,create_boundary(S),'T',0);
Ksto = calc_rigipc(S,X);

f=bodyload(S,[],'QN',1);

fsto = f.*X{4};


%% Calcul PCG
tic
qpc = cgs(Ksto,fsto,1e-6);
toc

q = unfreevector(S,qpc);
%%
plot(FENODEFIELD(random(q)),S)

%% Direct-SD
nbfonc=8;
fsize=16;
qradref = spectral_decomposition(qpc,'reference',qpc,'nbfoncmax',nbfonc);
qradref = normsort(qradref);
figure(4)
clf
m=length(qradref);
if ~exist('optionsplot')
    optionsplot={'surface'};
end
nl=ceil(sqrt(m));nc=ceil(m/nl);
for i=1:m
    fullsubplot(nl,nc,i)
    text(1.5,0.3,['U_{' num2str(i) '}'],'FontSize',fsize)
    plot(FENODEFIELD(qradref{i}),S,optionsplot{:})
    axis off
end


%% PU-GSD
G = GSDSOLVER('tol',1e-3,'pfixmax',3,'pfixtol',5e-2,...
    'display',true,'update',true,'nbfoncmax',nbfonc);
tic;qrad = solve(G,Ksto,fsto,[],'reference',qpc);toc
qrad = normsort(qrad);
%tic;qrad = funsolve(G,@(a,b) calc_rigiexpect(S,a,b),Ksto,fsto,[],'reference',qpc);toc;
figure(2)
clf
m=length(qrad);
if ~exist('optionsplot')
    optionsplot={'surface'};
end
nl=ceil(sqrt(m));nc=ceil(m/nl);
for i=1:m
    fullsubplot(nl,nc,i)
    text(1.5,0.3,['U_{' num2str(i) '}'],'FontSize',fsize)
    plot(FENODEFIELD(qrad{i}),S,optionsplot{:})
    axis off
end

%% A-GSD
G = GSDSOLVER('type','arnoldi','tol',1e-5,'pfixmax',3,'pfixtol',5e-2,...
    'display',true,'update',true,'nbfoncmax',nbfonc,'restart',4);
qrad = solve(G,Ksto,fsto,[],'reference',qpc);
qrad = normsort(qrad);
%tic;qrad = funsolve(G,@(a,b) calc_rigiexpect(S,a,b),Ksto,fsto,[],'reference',qpc);toc;
figure(2)
clf
m=length(qrad);
if ~exist('optionsplot')
    optionsplot={'surface'};
end
nl=ceil(sqrt(m));nc=ceil(m/nl);
for i=1:m
    fullsubplot(nl,nc,i)
    text(1.5,0.3,['U_{' num2str(i) '}'],'FontSize',fsize)
    plot(FENODEFIELD(qrad{i}),S,optionsplot{:})
    axis off
end


%% A-GSD
qrada = solve_spectral_arnoldi(Ksto,fsto,'display',...
    'nbfoncmax',nbfonc,'restart',4,'reference',qpc);
qrada = normsort(qrada);
figure(2)
clf
m=length(qrada);
if ~exist('optionsplot')
    optionsplot={'surface'};
end
nl=ceil(sqrt(m));nc=ceil(m/nl);
for i=1:m
    fullsubplot(nl,nc,i)
    text(1.5,0.3,['U_{' num2str(i) '}'],'FontSize',fsize)
    plot(FENODEFIELD(qrada{i}),S,optionsplot{:})
    axis off
end

%%
