function compshell3Dpout2(elemtype,ddlforce,color,ampl)
% function compshell3Dpout2(elemtype,ddlforce,color,ampl)
% elemtype : element type ('DKT', 'DKQ',  'DST', 'DSQ', 'COQ4') ('DKQ' by default)
% ddlforce : degree of freedom for loading ('FX', 'FY', 'FZ') ('FZ' by default)
% color : color for display ('b' by default)
% ampl : amplification for display (L/max(abs(q)/4) by default)

if nargin<1 || isempty(elemtype)
    elemtype = 'DKQ';
end
if nargin<2 || isempty(ddlforce)
    ddlforce = 'FZ';
end
if nargin<3 || isempty(color)
    color = 'b';
end

rb = 3;
rs = 10;
H = 1;
L = 1;
l = 1;

db = 0.01;
d = 0.01;
IY = pi*db^4/2;
IZ = IY;

mat=MATERIALS();
mat{1} = ELAS_BEAM('E',70e9,'NU',0.3,'S',d*l,'IZ',IZ,'IY',IY,'IX',IY+IZ,'RHO',7600);
mat{2} = ELAS_SHELL('E',70e9,'NU',0.3,'DIM3',d);

P1 = POINT([0,0,H]);
P2 = POINT([L,0,H]);
P3 = POINT([L,l,H]);
P4 = POINT([0,l,H]);
P5 = POINT([0,l,0]);
P6 = POINT([L,l,0]);

D1 = DROITE(P1,P2);
D2 = DROITE(P2,P3);
D3 = DROITE(P3,P4);
D4 = DROITE(P4,P1);

S1 = mesh(LIGNE(P5,P4),rb,mat{1});
S1 = convertelem(S1,'BEAM');
S3 = mesh(LIGNE(P6,P3),rb,mat{1});
S3 = convertelem(S3,'BEAM');
S2 = mesh(QUADRANGLE(P1,P2,P3,P4),2*rs,2*rs,mat{2});
S2 = convertelem(S2,elemtype,'fullintegration');
S = union(S1,S2,S3);

S = final(S);
S = addcl(S,[P6;P5],{'U','R'},0);

F = surfload(S,D1,ddlforce,1e5);

K = calc_rigi(S);

q = K\F;

if nargin<4 || isempty(ampl)
    ampl=L/max(abs(q))/4;
end

figure(2)
plot(S,'Color','k','FaceColor','k','FaceAlpha',0.1)
plot(S+ampl*q,'Color',color,'FaceColor',color,'FaceAlpha',0.3)
lighting gouraud
axis on
xlabel('x')
ylabel('y')
