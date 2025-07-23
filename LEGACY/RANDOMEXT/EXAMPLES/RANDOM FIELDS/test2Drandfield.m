r=10;
e=0.5;L=1;H=1;
r1=r;r2=r*(L-e)/e;

P1=POINT([0,0]);
P8=POINT([0,H-e]);
P2=POINT([e,H]);
P3=POINT([e,H-e]);
P4=POINT([L-e,H]);
P5=POINT([L-e,0]);
P6=POINT([L,H]);
P7=POINT([L,0]);


S1=mesh(DOMAIN(2,P1,P3),r1,r2);
S2=mesh(DOMAIN(2,P8,P2),r1,r1);
S3=mesh(DOMAIN(2,P3,P6),r2,r1);
S = union(S1,S2,S3);
S = concatgroupelem(S);

m=4;
E = RANDFIELD(RFMARGINAL(RVUNIFORM(0.5,1.5)),EXP2CORREL(e));
E = KL(E,m,S,'order',6);
E = spectral_decomposition(E,'nbfoncmax',10)

type = 'ther';

mat=MATERIALS();
if strcmp(type,'meca')    
    mat{1}=ELAS_ISOT('E',1,'NU',0.3,'DIM3',1,'RHO',1);
else
    mat{1}=FOUR_ISOT('k',FENODEFIELD(E{1}));
end
S=setmaterial(S,mat{1});

figure(5)
clf
for i=1:m
    subplot(ceil(sqrt(m)),ceil(m/ceil(sqrt(m))),i)
    plot(FENODEFIELD(E{i}),S)
    axis off
    axis image
    pause(0.01)
end


S=final(S);


D1  = STRAIGHTLINE(P1,P7);
D3  = STRAIGHTLINE(P6,VECTOR([1;0]));

if strcmp(type,'meca')
    S=addcl(S,D1,{'U'},[0;0]);
    f=surfload(S,D3,'FY',-1);
    K=calc_rigi(S);
    q=K\f;
    figure
    clf

    ampl=0.01; 

    s=calc_sigma(S,q,'smooth');
    plot(s,S+ampl*q,'compo','SMYY','edges')
else
    S=addcl(S,[],{'T'},[0]);
    f=bodyload(S,[],'QN',1);
    K=calc_rigi(S);
    q=K\f;
    figure
    clf

    plot(FENODEFIELD(unfreevector(S,q)),S)
end


