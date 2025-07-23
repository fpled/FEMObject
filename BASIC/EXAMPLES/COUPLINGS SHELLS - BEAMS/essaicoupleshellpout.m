
%MATERIAUX
Ipieu = 5.14e-3;
Spieu = 0.1233 ;
Epieu = 2.1e11;
Edalle = 4e10;
Stirant = Spieu;
mat=MATERIALS();
mat{1} = ELAS_BEAM('E',Epieu,'NU',0.3,...
    'S',Spieu,'IZ',Ipieu,'IY',Ipieu,'IX',2*Ipieu);
mat{2} = ELAS_SHELL('E',4e10,'NU',0.3,'DIM3',0.3);
mat{3} = ELAS_ISOT('E',Epieu,'NU',0.3,'S',Stirant);

nbtroncons = 2;
l = 10;
L = l*nbtroncons ;
h = 10;
pasL = l;
Ppieu1 = POINT([0,0,0;pasL,0,0;pasL,l,0;0,l,0]);

profpieu = [-l,-l,-l,-l];

Ppieu2 = Ppieu1;
for i=1:4
    Ppieu2(i) = Ppieu1(i)+VECTOR([0;0;profpieu(i)]);
end

rpieu = 7;
Mpieu = MODEL('TRID');
for i=1:4
    Mpieu = union(Mpieu,mesh(LINE(Ppieu1(i),Ppieu2(i)),rpieu,mat{1}));
end
Mpieu = concatgroupelem(Mpieu);
Mpieu = convertelem(Mpieu,'BEAM');
r1dalle = 5*3;
r2dalle = 5*3;
Pdalle = POINT([0,0,0;pasL,0,0;pasL,l,0;0,l,0]);
Mdalle = mesh(QUADRANGLE(Pdalle(1),Pdalle(2),...
    Pdalle(3),Pdalle(4)),r1dalle,r2dalle,mat{2});
Mdalle = convertelem(Mdalle,'DKQ'); %% on peut aussi utiliser COQ4

M = union(Mpieu,Mdalle);
figure(1)
clf
plot(M,'FaceColor','w','node')

V = VECTOR([pasL;0;0]);
Mquai = M;
Pbaspieu = Ppieu2;
for k=1:nbtroncons-1
    M = M+V;
    Ppieu2 = Ppieu2 + V ;
    Mquai = union(Mquai,M);
    Pbaspieu = [Pbaspieu,Ppieu2];
end

figure(2)
clf
plot(Mquai,'FaceColor','w')

Mquai = final(Mquai);
Mquai = addcl(Mquai,Pbaspieu,{'U','R'});

D = STRAIGHTLINE(POINT([0,0,0]),POINT([L,0,0]));
F = surfload(Mquai,D,'FX',-100000);
F = F + surfload(Mquai,D,'FY',-200000);

K = calc_rigi(Mquai);

q = K\F;

figure(3)
clf
ampl = 10;
plot(Mquai,'Color','w','FaceColor','w','FaceAlpha',0.3);
plot(Mquai+ampl*q,'Color','b','FaceColor','b','FaceAlpha',0.3);
