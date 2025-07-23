function shell3D(test,elemtype,ddlforce,color,ampl)
% function shell3D(test,elemtype,ddlforce,color,ampl)
% test = 1 : 1 plaque, chargement lineique de flexion aux bords de la plaque
% test = 2 : assemblage 2 plaques, chargement lineique a la jonction des 2 plaques (by default)
% test = 3 : 1 plaque, chargement surfacique de flexion sur la plaque
% test = 4 : assemblage 2 plaques, chargement nodal a la jonction des 2 plaques
% elemtype : element type ('DKT', 'DKQ', 'DST', 'DSQ', 'COQ4') ('DKQ' by default)
% ddlforce : degree of freedom for loading ('FX', 'FY', 'FZ') ('FZ' by default)
% color : color for display ('b' by default)
% ampl : amplification for display (L/max(abs(q)) by default)

if nargin<1 || isempty(test)
    test = 2;
end
if nargin<2 || isempty(elemtype)
    elemtype = 'DKQ';
end
if nargin<3 || isempty(color)
    color = 'b';
end
if nargin<4 || isempty(ddlforce)
    ddlforce = 'FZ';
end

r = 20;
H = 1;
L = 1;
l = 1;

dim3 = 0.01;
mat = MATERIALS();
mat{1} = ELAS_SHELL('E',70e9,'NU',0.3,'DIM3',dim3,'RHO',7600);

P1 = POINT([0,0,0]);
P2 = POINT([L,0,0]);
P3 = POINT([L,l,0]);
P4 = POINT([0,l,0]);
P5 = POINT([L,0,H]);
P6 = POINT([L,l,H]);

D1 = STRAIGHTLINE(P1,P4);
D2 = STRAIGHTLINE(P2,P3);
D3 = STRAIGHTLINE(P5,P6);

S = mesh(QUADRANGLE(P1,P2,P3,P4),r,r,mat{1});
S = convertelem(S,elemtype);

if test==2 || test==4
    S2 = mesh(QUADRANGLE(P2,P3,P6,P5),r,r,mat{1});
    S2 = convertelem(S2,elemtype);
    S = union(S,S2);
end

S = final(S);

S = addcl(S,D1,{'U','R'},0);

if test==1
    F = surfload(S,D2,ddlforce,1e5);
elseif test==2
    F = surfload(S,D3,ddlforce,1e5);    
elseif test==3
    F = bodyload(S,[],ddlforce,1e5);
elseif test==4
    F = nodalload(S,D3,ddlforce,1e5);
end
K = calc_rigi(S);
M = calc_mass(S);

q = K\F;

if nargin<5 || isempty(ampl)
    ampl = L/max(abs(q));
end

figure(1)
clf
plot(S,'Color','k','FaceColor','k','FaceAlpha',0.1)
plot(S+ampl*q,'Color',color,'FaceColor',color,'FaceAlpha',0.3)
camlight left
lighting gouraud

figure(2)
clf
plot(S)
hold on
[hD,legD] = plotbcond(S);
[hN,legN] = vectorplot(S,'F',F,5*ampl,'r','LineWidth',1);
legend([hD,hN],[legD,legN])

P23 = POINT([L,l/2,0]);
P56 = POINT([L,l/2,H]);

q1z = eval_sol(S,q,P1,'UZ')
q2z = eval_sol(S,q,P2,'UZ')
q3z = eval_sol(S,q,P3,'UZ')
q4z = eval_sol(S,q,P4,'UZ')
q23z  = eval_sol(S,q,P23,'UZ')
if test==2 || test==4
    q56ux = eval_sol(S,q,P56,'UX')
    q56uy = eval_sol(S,q,P56,'UY')
    q56uz = eval_sol(S,q,P56,'UZ')
    q56rx = eval_sol(S,q,P56,'RX')
    q56ry = eval_sol(S,q,P56,'RY')
    q56rz = eval_sol(S,q,P56,'RZ')
end
