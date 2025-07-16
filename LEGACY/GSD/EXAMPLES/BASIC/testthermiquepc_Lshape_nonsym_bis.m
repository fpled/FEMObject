%% modele stochastique
r=13;
p=4;

P = POINT([ 0,0 ; 1,0 ; 0,2 ; 1,2 ; 1,1 ; 2,1 ; 2,2 ]);
S1 = mesh(DOMAIN(2,P(1),P(4)),r,r*2);
S2 = mesh(DOMAIN(2,P(5),P(7)),r,r);
S = union(S1,S2);
S=concatgroupelem(S);
secdet=1;

RV=RANDVARS();

RV{1}=RVUNIFORM(0.1,2);
RV{2}=RVUNIFORM(-3,3);
RV{3}=RVUNIFORM(0,6);

if ~secdet
    RV{4}=RVNORMAL(0.5,0.2);
    RV{5}=RVNORMAL(0,0.2);
end

X = PCMODEL(RV,'order',p,RANDPOLYS(RV));
PC = getPC(X);
mat=FOUR_ISOT('k',X{1},'b',[-X{2};-X{3}]);
S = setmaterial(S,mat);
S=final(S);
L1 = LIGNE(P(6),P(7));
L2 = LIGNE(P(1),P(3));
L3 = LIGNE(P(3),P(7));
L4 = LIGNE(P(6),P(5));
L5 = LIGNE(P(2),P(5));
L6 = LIGNE(P(1),P(2));

% conditions aux limites
S=addcl(S,L1,'T',0);
S=addcl(S,L2,'T',0);
S=addcl(S,L3,'T',0);
S=addcl(S,L4,'T',0);
S=addcl(S,L5,'T',0);


f1=bodyload(S,[],'QN',1);
f2 = surfload(S,L6,'QN',1);
if secdet
    fsto = f1+f2;    
else
    fsto = f1*X{end-1}+f2.*X{end};
end


%mat=FOUR_ISOT('k',1,'b',0*[1;1]);
%S=setmaterial(S,mat);
%K1 = calc_rigi(S);
%mat=FOUR_ISOT('k',0,'b',[-1;0]);
%S=setmaterial(S,mat);
%K2 = calc_rigi(S);
%mat=FOUR_ISOT('k',0,'b',[0;-1]);
%S=setmaterial(S,mat);
%K3 = calc_rigi(S);

%Ksto = K1.*X{1}+K2.*X{2}+K3.*X{3};
Ksto = calc_rigipc(S,PC);


%% resolution CGS
tic
qpc = cgs(Ksto,fsto,1e-8);
toc

q = unfreevector(S,qpc);
plot(FENODEFIELD(random(q)),S)


nbfonc=8;
fsize=16;
pfixmax=4;
pfixtol=5e-2;

[qradref,result_direct] = spectral_decomposition(qpc,'reference',qpc,'nbfoncmax',nbfonc,'testerror');
qradref = normsort(qradref);

figure(21)
clf
nl=3;nc=3;
for k=1:nbfonc
    subplot(nl,nc,k)

    plot(FENODEFIELD(unfreevector(S,qradref{k})),S)
    axis off
    text(1.5,0.3,['U_{' num2str(k) '}'],'FontSize',fsize)
end

%%
GSD = GSDSOLVER('tol',1e-7,'nbfoncmax',20,...
    'display',true,'update',true,'errorindicator','rayleigh','type',...
    'power','restart',0,'orthocrit',1e-12,'nbfoncmaxsimul',12);
GSD = setparam(GSD,'errorindicator','reference');

figure(100)
clf
% PGSD 
GSD = setparam(GSD,'update',false);
GSD = setparam(GSD,'finalupdate',true);
GSD = setparam(GSD,'orthoduringpfix',0);
[qradpugsd,resultpgsd]=solve(GSD,Ksto,fsto,[],'reference',qpc);
semilogy(resultpgsd.error,'b-*')
hold on
% PUGSD 
GSD = setparam(GSD,'update',true);
GSD = setparam(GSD,'orthoduringpfix',0);
GSD = setparam(GSD,'finalupdate',false);
[qradpugsd,resultpugsd]=solve(GSD,Ksto,fsto,[],'reference',qpc);
semilogy(resultpugsd.error,'k->')
hold on
% PGSD avec reortho interne
GSD = setparam(GSD,'orthoduringpfix',1);
GSD = setparam(GSD,'update',false);
GSD = setparam(GSD,'finalupdate',true);
[qradpugsd,resultpgsdort]=solve(GSD,Ksto,fsto,[],'reference',qpc);
semilogy(resultpgsdort.error,'r-<')
% PUGSD avec reortho interne
GSD = setparam(GSD,'orthoduringpfix',1);
GSD = setparam(GSD,'update',true);
GSD = setparam(GSD,'finalupdate',false);
[qradpugsd,resultpugsdort]=solve(GSD,Ksto,fsto,[],'reference',qpc);
semilogy(resultpugsdort.error,'m-v')


GSD = setparam(GSD,'type','arnoldi');
GSD = setparam(GSD,'orthocrit',1e-12); 
GSD = setparam(GSD,'nbfoncmaxsimul',10); 
GSD = setparam(GSD,'restart',1); 
[qradpugsd,resultagsd]=solve(GSD,Ksto,fsto,[],'reference',qpc);

semilogy(resultagsd.error,'k--d')

% PUGSD - SUBSPACE ITERATIONS
GSD = setparam(GSD,'type','power');
GSD = setparam(GSD,'update',true);
GSD = setparam(GSD,'orthoduringpfix',0);
GSD = setparam(GSD,'finalupdate',false);
GSD = setparam(GSD,'subspaceiteration',2);
GSD = setparam(GSD,'nbfoncmax',9);
[qradpugsdsub,resultpugsdsub]=solve(GSD,Ksto,fsto,[],'reference',qpc);
semilogy(resultpugsdsub.error,'y-o')
hold on

% PUGSD - SUBSPACE ITERATIONS - RESIDUAL
GSD = setparam(GSD,'type','power');
GSD = setparam(GSD,'update',true);
GSD = setparam(GSD,'orthoduringpfix',0);
GSD = setparam(GSD,'finalupdate',false);
GSD = setparam(GSD,'subspaceiteration',3);
GSD = setparam(GSD,'nbfoncmax',9);
Kmsym = expect(Ksto)/2+expect(Ksto)'/2;
KmKsto = (Kmsym\Ksto);
[qradpugsdres,resultpugsdres]=solve(GSD,KmKsto'*Ksto,KmKsto'*fsto,[],'reference',qpc);
semilogy(resultpugsdres.error,'g-s')
hold on



legend('PGSD avec update final','PUGSD','PGSD avec reortho k et final update','PUGSD avec reortho k','arnoldi','PUGSD - subspaceiterations','PUGSD - residual')
