r=100;
L = 1;
P=POINT([0;L]);
M = mesh(LINE(P(1),P(2)),r);

RV = RANDVARS();
RV{1}=RVUNIFORM(0.5,1.5);
RV{2}=RVNORMAL(1,0.2);
RV{3}=RVBETA(3,3,1,2);

X = PCMODEL(RV,'order',3,RANDPOLYS(RV));

mat=MATERIALS();
mat{1}=ELAS_ISOT('E',1,'S',1,'NU',0.3);

M=setmaterial(M,mat{1});
M=final(M);

M=addcl(M,P(1),'U',0);

clear K
K{1} = calc_rigi(M);

clear f
f{1} = nodalload(M,P(2),'FX',1);
f{2} = bodyload(M,[],'FX',1);
Ksto = K{1}.*X{1};
fsto = f{1}.*X{2} + f{2}*X{3};

ampl = 1 ;

Kmean = mean(Ksto) ;
fmean = mean(fsto) ;
u0=Kmean\fmean;
plot(M+u0*ampl,'EdgeColor','r');


RVsample = random(RV);
Ksample = randomeval(Ksto,RVsample);
fsample = randomeval(fsto,RVsample);
usample=Ksample\fsample;
plot(M+usample*ampl,'EdgeColor','b');


upc = pcg(Ksto,fsto,1e-10);

figure(1)
plot(M+mean(upc)*ampl,'node')

figure(2)
s=calc_sigma(M,mean(upc));
plot(s,M+mean(upc)*ampl,'compo','SMXX','node')

[urad,result2]=solve_spectral(Ksto,fsto,'nbfoncmax',9,'reference',upc);
[urad2,result2]=solve_spectral_arnoldi(Ksto,fsto,'nbreact',0,'nbfoncmaxsimul',7,'reference',upc);
