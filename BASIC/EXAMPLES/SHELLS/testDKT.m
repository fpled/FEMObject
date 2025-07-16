function testDKT(elemtype,ddlforce,color,ampl)
% function testDKT(elemtype,ddlforce,color,ampl)
% elemtype : element type ('DKT', 'DKQ', 'DST', 'DSQ', 'COQ4', 'TRI3', 'QUA4') ('DKQ' by default)
% ddlforce : degree of freedom for loading ('FX', 'FY', 'FZ') ('FZ' by default)
% color : color for display ('b' by default)
% ampl : amplification for display (L/max(abs(q)) by default)

if nargin<1 || isempty(elemtype)
    elemtype = 'DKQ';
end
if nargin<2 || isempty(ddlforce)
    ddlforce = 'FZ';
end
if nargin<3 || isempty(color)
    color = 'b';
end

r = 10;
H = 1;
L = 1;
l = 1/3;

dim3 = 0.01;
mat = MATERIALS();
if strcmp(elemtype,'TRI3') || strcmp(elemtype,'QUA4')
    mat{1} = ELAS_ISOT('E',70e9,'NU',0.3,'DIM3',dim3,'RHO',7600);
    P1=POINT([0,0]);
    P2=POINT([L,0]);
    P3=POINT([L,l]);
    P4=POINT([0,l]);
    P5=POINT([L,l/2]);
else
    mat{1} = ELAS_SHELL('E',70e9,'NU',0.3,'DIM3',dim3,'RHO',7600);
    P1=POINT([0,0,0]);
    P2=POINT([L,0,0]);
    P3=POINT([L,l,0]);
    P4=POINT([0,l,0]);
    P5=POINT([L,l/2,0]);
end

D1 = DROITE(P1,P4);
D2 = DROITE(P2,P3);

S = mesh(QUADRANGLE(P1,P2,P3,P4),r,r,mat{1});
S = convertelem(S,elemtype);
S = final(S);

S = addcl(S,D1,{'U','R'},0);

% I = l^3/12;
% alpha = 1e5*l*H/I;
% % fun = inline('-(x(:,2)-l/2)*alpha','x','l','alpha');
% fun = @(x,l,alpha) -(x(:,2)-l/2)*alpha;
% F = surfload(S,D2,ddlforce,fun,l,alpha);
F = surfload(S,D2,ddlforce,1);
% F = bodyload(S,[],ddlforce,1e5);
% F = nodalload(S,D2,ddlforce,1e5);
K = calc_rigi(S);
M = calc_mass(S);

q = K\F;

if nargin<4 || isempty(ampl)
    if strcmp(elemtype,'TRI3') || strcmp(elemtype,'QUA4')
        ampl = L/max(abs(q))/10;
    else
        q = unfreevector(S,q);
        Q = q(findddl(S,DDL(DDLVECT('U',S.syscoord,'TRANS'))),:);
        ampl = L/max(abs(Q))/10;
    end
end

m = 1:10;

figure(length(m)+1)
clf
plot(S,'Color','k','FaceColor','k','FaceAlpha',0.1)
plot(S+ampl*q,'Color',color,'FaceColor',color,'FaceAlpha',0.3)
if strcmp(ddlforce,'FZ')
    view(3)
else
    view(2)  
end
camlight left
lighting gouraud

q2x = eval_sol(S,q,P2,'UX')
q2y = eval_sol(S,q,P2,'UY')
if ~(strcmp(elemtype,'TRI3') || strcmp(elemtype,'QUA4'))
    q2z = eval_sol(S,q,P2,'UZ')
end

[V,D] = calc_mode(K,M,m);

for i=m
    figure(i)
    clf
    if strcmp(elemtype,'TRI3') || strcmp(elemtype,'QUA4')
        ampl = L/max(abs(V(:,i)))/10;
    else
        vi = unfreevector(S,V(:,i));
        Vi = vi(findddl(S,DDL(DDLVECT('U',S.syscoord,'TRANS'))),:);
        ampl = L/max(abs(Vi))/10;
    end
    plot(S,'Color','k','FaceColor','k','FaceAlpha',0.1)
    plot(S+ampl*V(:,i),'Color',color,'FaceColor',color,'FaceAlpha',0.3)
    camlight left
    lighting gouraud
end
