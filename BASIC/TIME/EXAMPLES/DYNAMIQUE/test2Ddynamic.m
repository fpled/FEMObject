

% D = mesh(DOMAIN(2,P(1),P(3)),2*n,2*n);
% DD = mesh(DOMAIN(2,P(5),P(7)),n,n);
% D = union(D,DD);
% S = MODEL('PLANE');
% S = addnode(S,D.node);
% S = addelem(S,'QUA4',D,'option','CONT','mat',MAT);

MAT = ELAS_ISOT('E',1,'NU',0.,'RHO',1,'DIM3',1);
probtype = 1;
switch probtype
    case 1
        D = DOMAIN(2,[0,0],[2,.3]);
        S = gmsh(D,.03);
        
        S = setmaterial(S,MAT);
    case 2
        S = cast2matlab_model('PLANE','C:\PROGRAMMES\CASTEM\EXAMPLES\mesh2Ddyna_model.txt',MAT);
        % S = cast2matlab_model('PLANE','E:\PROGRAMMES\CASTEM\EXAMPLES\mesh2Ddyna_model_gros_1.txt',MAT);
end
xnode = getcoord(S.node);
xmax = max(xnode(:,1));
xmin = min(xnode(:,1));
ymax = max(xnode(:,2));
ymin = min(xnode(:,2));

P = POINT([xmin,ymin;xmax,ymin;xmax,ymax;xmin,ymax]);

D1 = LINE(P(1),P(4));
D2 = LINE(P(2),P(3));

% S = convertelem(S,'TRI3');
S = final(S);
S = addcl(S,D1,'U',[0;0]);
M = calc_mass(S);
K = calc_rigi(S);
f = surfload(S,D2,'FX',-1);

if ~exist('t0','var')
    t0 = 0;
end
if ~exist('t1','var')
    t1 = 2;
end
if ~exist('nt','var')
    nt = 50;
end

T = TIMEMODEL(t0,t1,nt);

if ~exist('tc','var')
    tc = get(T,'t1')/6;
end


if ~exist('solv','var')
    solv = 2;
end
if ~exist('psolv','var')
    psolv = 1;
end

switch solv
    case 1
        N = NEWMARKSOLVER(T,'alpha',0.05);
    case 2
        N = DGTIMESOLVER(T,psolv);
    otherwise
        error('pas defini')
end
N = setparam(N,'display',true);

if ~exist('loadcase','var')
    loadcase = 1;
end

switch loadcase
    case 1
        ft = f*rampe(N,t0,tc);
    case 2
        ft = f*dirac(N,t0,tc);
    case 3
        ft = f*one(N);
end

%% resolution
[ut,result,vt,at] = ddsolve(N,ft,M,K);

em = calc_epsilon(S,ut);
sm = calc_sigma(S,ut);

%% affichage solution
ut = setevolparam(ut,'colormap',jet,'colorbar',true);
vt = setevolparam(vt,'colormap',jet,'colorbar',true);
at = setevolparam(at,'colormap',jet,'colorbar',true);

figure(1)
clf
evol_sol(ut,S,'displ',1,'rescale',true)

figure(2)
clf
evol_sol(vt,S,'displ',1,'rescale',true)

figure(3)
clf
evol_sol(at,S,'displ',1,'rescale',true)

figure(4)
clf
evol_sol(ut,S,'sigma',1,'rescale',true)

smmax = max(max(sm(1)));
smmin = min(min(sm(1)));
N = setevolparam(N,'plotstep',1,'caxis',[smmin,smmax],...
    'pausetime',1/nt,'colormap',jet,'colorbar',true);

figure(5)
clf
evol(N,sm,S,'compo','SMXX')

N = setevolparam(N,'view',3,'axison',true,'axis',[xmin,xmax,ymin,ymax,smmin,smmax]);

figure(6)
clf
evol(N,sm,S,'compo','SMXX')

figure(7)
clf
optionsplot = {'k'};
hold on
se = getvalue(sm);
se = se{1};
dx = xmax-xmin;
dy = ymin-ymax;
[NT,repnode] = getnodenextto(S.node,POINT([dx/2,dy/2]),'local');
plot(N,double(se(1,:,repnode)),optionsplot{:});
