nbtroncons = 37;
nbblocelem = 1;
rpieu = 6;
r1dalle = 2*1;
r2dalle = 7*1;
rtirant = 1;

%MATERIAUX
Ipieu = 5.14e-3;
Spieu = 0.1233 ;
Epieu = 2.1e11;
Edalle = 4e10;
Stirant = 5.675e-3;

mat=MATERIALS();
mat{1} = ELAS_BEAM('E',Epieu,'NU',0.3,...
    'S',Spieu,'IZ',Ipieu,'IY',Ipieu,'IX',2*Ipieu);
mat{2} = ELAS_SHELL('E',4e10,'NU',0.3,'DIM3',0.3);
mat{3} = ELAS_ISOT('E',1,'NU',0.3,'S',Stirant);

L = 6.68*nbtroncons; 
l = 42; 
ltirant = 20;
proftirant = -1.2;
profpieu = [-21,-19,-17,-15,-13,-11,-9,-7];

pasL = L/nbtroncons/2;
pasl = l/7;

Ppieu1 = POINT([0,0,0;pasl,pasL,0;2*pasl,0,0;3*pasl,pasL,0;...
    4*pasl,0,0;5*pasl,pasL,0;6*pasl,0,0;7*pasl,pasL,0]);
Ppieu2 = Ppieu1;
Ptirant1 = POINT([l,pasL,0]);
Ptirant2 = POINT([l+ltirant,pasL,proftirant]);
for i=1:8
    Ppieu2(i) = Ppieu1(i)+VECTOR([0;0;profpieu(i)]);
end

Mpieu = MODEL('TRID');
for i=1:8
    Mpieu = union(Mpieu,mesh(LINE(Ppieu1(i),Ppieu2(i)),rpieu,mat{1}));    
end
Mpieu = concatgroupelem(Mpieu);
Mpieu = convertelem(Mpieu,'BEAM');
Mtirant = mesh(LINE(Ptirant1,Ptirant2),rtirant,mat{3});
Mtirant = convertelem(Mtirant,'BARR');

Pdalle = POINT([0,0,0;0,2*pasL,0;l,2*pasL,0;l,0,0]);
Mdalle = mesh(QUADRANGLE(Pdalle(1),Pdalle(2),...
    Pdalle(3),Pdalle(4)),r1dalle,r2dalle,mat{2});
Mdalle = convertelem(Mdalle,'DKQ');

M = union(Mpieu,Mdalle);
M = union(M,Mtirant);

figure(1)
clf
plot(M,'FaceColor','w','node')

if nbblocelem>1
    V = VECTOR([0;2*pasL;0]);
    Mquai = M;
    Pbastirant = Ptirant2;
    Pbaspieu = Ppieu2;
    for k=1:nbblocelem-1
        M = M+V;
        Ppieu2 = Ppieu2 + V ;
        Ptirant2 = Ptirant2 + V ;
        Mquai = union(Mquai,M);
        Pbaspieu = [Pbaspieu,Ppieu2];
        Pbastirant = [Pbastirant,Ptirant2];
    end
    Ppieu2 = Pbaspieu;
    Ptirant2 = Pbastirant;
    M = Mquai;

    figure(2)
    clf
    plot(Mquai,'group')
    pause(0.001)
    M = concatgroupelem(M);
end

V = VECTOR([0;nbblocelem*2*pasL;0]);
Mquai = M;
Pbastirant = Ptirant2;
Pbaspieu = Ppieu2;
for k=1:nbtroncons/nbblocelem-1
    M = M+V;
    Ppieu2 = Ppieu2 + V ;
    Ptirant2 = Ptirant2 + V ;
    Mquai = union(Mquai,M);
    Pbaspieu = [Pbaspieu,Ppieu2];
    Pbastirant = [Pbastirant,Ptirant2];
end

figure(2)
clf
plot(Mquai,'FaceColor','w')
Mquai = concatgroupelem(Mquai);
Mquai = final(Mquai);

Mquai = addcl(Mquai,Pbaspieu,{'U','R'});
Mquai = addcl(Mquai,Pbastirant,{'U'});

D = STRAIGHTLINE(POINT([0,0,0]),POINT([0,L,0]));
F = surfload(Mquai,D,'FX',-100000);

RV = RANDVARS(RVLOGNORMAL(Epieu,Epieu*0.3,'stat'));
X = PCMODEL(RV,'order',5);

Kd = calc_rigi(Mquai,'selgroup',1);
Kp = calc_rigi(Mquai,'selgroup',2);
Kt = calc_rigi(Mquai,'selgroup',3);

K = Kd+Kp+Kt*X{1};

qpc = pcg(K,F,1e-6);

q = mean(qpc);
figure(3)
clf
ampl = 1000;
plot(Mquai,'Color','w');
plot(Mquai+ampl*q,'Color','r','FaceColor','r','FaceAlpha',0.5);

qbout = eval_sol(Mquai,qpc,POINT([0,0,0]),'UX')

figure(4)
clf
deltaeau = -4;
ploteausol = 1;

if ploteausol==1
    eau = QUADRANGLE(POINT([-l/5,-L/10,deltaeau]),POINT([l,-L/10,deltaeau]),POINT([l,L+L/10,deltaeau]),POINT([-l/5,L+L/10,deltaeau]));
    plot(eau,'FaceColor','b','EdgeColor','none','FaceAlpha',0.3)
    sol = QUADRANGLE(POINT([0,-L/10,-21]),POINT([l,-L/10,-7]),POINT([l,L+L/10,-7]),POINT([0,L+L/10,-21]));
    plot(sol,'FaceColor','r','EdgeColor','none','FaceAlpha',0.3)
    camlight headlight;
    lighting gouraud
end

plot(Mquai,'FaceColor','m','EdgeColor','m','FaceAlpha',0.3)
camlight headlight;
lighting gouraud
