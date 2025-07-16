
L = 60.;
H = 30.;
H1= 10.;
H2 = 20;
B = 5 ;

M = cast2matlab_model('PLAN','E:\PROGRAMMES\CASTEM\EXAMPLES\foundationfine_model.txt');

RV = RANDVARS();
RV{1}=RVBETA(3,3,30,70);
RV{2}=RVBETA(3/2,3/2,30,70);
RV{3}=RVUNIFORM(30,70);
RV{4}=RVLOGNORMAL(0.2,0.2*0.25,'stat');
RV{5}=RVGAMMA(0.2,0.2*0.25,'stat',0.1);


X = PCMODEL(RV,'order',3,RANDPOLYS(RV));


mat=MATERIALS();
mat{1}=ELAS_ISOT('E',X{1},'NU',0.3,'DIM3',1);
mat{2}=ELAS_ISOT('E',X{2},'NU',0.3,'DIM3',1);
mat{3}=ELAS_ISOT('E',X{3},'NU',0.3,'DIM3',1);

M=setmaterial(M,mat{1},1);
M=setmaterial(M,mat{2},2);
M=setmaterial(M,mat{3},3);
M=setoption(M,'CONT');


M=final(M);
plot(M)

M=addcl(M,DROITE(POINT([0,-H]),POINT([L,-H])),'U',0);
M=addcl(M,DROITE(POINT([0,0]),POINT([0,-H])),'UX',0);
M=addcl(M,DROITE(POINT([L,0]),POINT([L,-H])),'UX',0);

clear K

Ksto = calc_rigipc(M,X);

L1 = LIGNE(POINT([0,0]),POINT([B,0]));
L2 = LIGNE(POINT([B,0]),POINT([L,0]));
clear f
f{1} = surfload(M,L1,'FY',-1);
f{2} = surfload(M,L1,'FY',-1);


fsto = f{1}.*X{4}+f{2}.*X{5};

ampl = 20 ;

Kmean = mean(Ksto) ;
fmean = mean(fsto) ;
u0=Kmean\fmean;
plot(M+u0*ampl,'EdgeColor','r');

RVsample = random(RV);
Ksample = randomeval(Ksto,[RVsample{:}],RV);
fsample = f{1}*RVsample{4}+f{2}*RVsample{5};
usample=Ksample\fsample;
plot(M+usample*ampl,'EdgeColor','b');

tic
fprintf('Resolution pcg')
upc = pcg(Ksto,fsto,1e-7);
%upc = solve(Ksto,fsto);
fprintf('Temps pcg  ')
toc
figure
plot(M+mean(upc)*ampl,'EdgeColor','g')

%[urad,radU,radl,radalpha,result]=solve_radial_simple_gc(Ksto,fsto,paramresol.radial,double(upc));

[urad,result2]=solve_spectral(Ksto,fsto,'nbfoncmax',6,'pfixmax',6,'pfixtol',1e-3);

[urad2,result2]=solve_spectral_arnoldi(Ksto,fsto,'nbreact',0,'orthocrit',1e-12,'nbfoncmaxsimul',6,'reference',upc);

for i=1:length(urad)
    epsrad{i} = calc_epsilon(M,urad{i});   
    figure(i)
    clf
    plot(epsrad{i}(2),M);
end


plot(M+mean(urad)*ampl,'EdgeColor','m')

