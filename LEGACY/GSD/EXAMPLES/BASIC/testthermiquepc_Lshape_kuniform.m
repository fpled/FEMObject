r=10;
P = POINT([ 0,0 ; 1,0 ; 0,2 ; 1,2 ; 1,1 ; 2,1 ; 2,2 ]);
S1 = mesh(DOMAIN(2,P(1),P(4)),r,r*2);
S2 = mesh(DOMAIN(2,P(5),P(7)),r,r);
S = union(S1,S2);

RV=RANDVARS();

secdet = 0;

e=0.4
k=RVUNIFORM(1-e,1+e)

if secdet
    RV = RANDVARS(RANDVARS(k));
else
    RV = RANDVARS(RANDVARS(k),RVNORMAL(0.5,0.2),RVNORMAL(0,0.2));
end
p=8;
X = PCMODEL(RV,'order',p,RANDPOLYS(RV));


L1 = LIGNE(P(6),P(7));
L2 = LIGNE(P(1),P(3));
L3 = LIGNE(P(3),P(7));
L4 = LIGNE(P(6),P(5));
L5 = LIGNE(P(2),P(5));
L6 = LIGNE(P(1),P(2));

mat=FOUR_ISOT('k',1);
S = setmaterial(S,mat);
S=final(S);
S=addcl(S,L1,'T',0);
S=addcl(S,L2,'T',0);
S=addcl(S,L3,'T',0);
S=addcl(S,L4,'T',0);
S=addcl(S,L5,'T',0);

mat=FOUR_ISOT('k',1);
S=setmaterial(S,mat);
K = calc_rigi(S);
Ksto = K*X{1};



f1=bodyload(S,[],'QN',1);
f2 = surfload(S,L6,'QN',1);
if secdet
    fsto = f1+f2;    
else
    fsto = f1*X{end-1}+f2.*X{end};
end

disp('Resolution pcg')
tic
if israndom(Ksto)
    qpc = pcg(Ksto,fsto,1e-14);
else
    qpc = Ksto\fsto;    
end
toc
time_pcg=toc;


figure(200)
clf
plot(S,'FaceColor',[0.6,0.9,0.9])
%plot(S)
colu = 'b' ;
colg = 'r' ;
linwid = 3;
plot(L1,'LineWidth',linwid,'EdgeColor',colu)
plot(L2,'LineWidth',linwid,'EdgeColor',colu)
plot(L3,'LineWidth',linwid,'EdgeColor',colu)
plot(L4,'LineWidth',linwid,'EdgeColor',colu)
plot(L5,'LineWidth',linwid,'EdgeColor',colu)
plot(L6,'LineWidth',linwid,'EdgeColor',colg)

annotation('textarrow',[0.7 0.65],[0.35 0.5],'String','\partial_1\Omega','FontSize',18)
annotation('textarrow',[0.55 0.45],[0.2 0.12],'String','\partial_2\Omega','FontSize',18)

%saveas(gcf,'E:\REDACTION\matlab_results\Lshape_mesh_boundaryconditions.eps','epsc2')
%saveas(gcf,'E:\REDACTION\matlab_results\Lshape_mesh_boundaryconditions_sanscol.eps','epsc2')


figure(201)
clf
%plot(S,'FaceColor',[0.6,0.9,0.9])
plot(S)
repnode=330;
plot(S.node(repnode),'bx','MarkerSize',12)
plottext(S.node(repnode),'P_1','Color','b','FontSize',18)
repnode=640;
plot(S.node(repnode),'bx','MarkerSize',12)
plottext(S.node(repnode),'P_2','FontSize',18,'Color','b')
repnode=1097;
plot(S.node(repnode),'bx','MarkerSize',12)
plottext(S.node(repnode),'P_3','FontSize',18,'Color','b')


%saveas(gcf,'E:\REDACTION\matlab_results\Lshape_mesh_boundaryconditions.eps','epsc2')
%saveas(gcf,'E:\REDACTION\matlab_results\Lshape_mesh_withpoints.eps','epsc2')

figure(202)
clf
plot(S,'FaceColor',[0.6,0.9,0.9])
%plot(S)
colu = 'b' ;
colg = 'r' ;
linwid = 3;
plot(L1,'LineWidth',linwid,'EdgeColor',colu)
plot(L2,'LineWidth',linwid,'EdgeColor',colu)
plot(L3,'LineWidth',linwid,'EdgeColor',colu)
plot(L4,'LineWidth',linwid,'EdgeColor',colu)
plot(L5,'LineWidth',linwid,'EdgeColor',colu)
plot(L6,'LineWidth',linwid,'EdgeColor',colg)

annotation('textarrow',[0.7 0.65],[0.35 0.5],'String','\partial_1\Omega','FontSize',18)
annotation('textarrow',[0.55 0.45],[0.2 0.12],'String','\partial_2\Omega','FontSize',18)

repnode=330;
plot(S.node(repnode),'Marker','kx','MarkerSize',12)
plottext(S.node(repnode),'P_1','Color','k','FontSize',18)
repnode=640;
plot(S.node(repnode),'Marker','kx','MarkerSize',12)
plottext(S.node(repnode),'P_2','FontSize',18,'Color','k')
repnode=1097;
plot(S.node(repnode),'Marker','kx','MarkerSize',12)
plottext(S.node(repnode),'P_3','FontSize',18,'Color','k')

%saveas(gcf,'E:\REDACTION\matlab_results\Lshape_mesh_boundaryconditions_withpoints.eps','epsc2')
%saveas(gcf,'E:\REDACTION\matlab_results\Lshape_mesh_boundaryconditions_sanscol.eps','epsc2')

%%
GSD = GSDSOLVER('tol',1e-4,'nbfoncmax',6,...
    'display',true,'update',true,'errorindicator','residual','type','arnoldi','restart',3);
[qrad,result]=solve(GSD,Ksto,fsto);

GSD = GSDSOLVER('tol',1e-4,'nbfoncmax',6,...
    'display',true,'update',true,'errorindicator','residual','type','power','restart',3);
[qrad,result]=solve(GSD,Ksto,fsto);

GSD = GSDSOLVER('tol',1e-4,'nbfoncmax',2,...
    'display',true,'update',true,'errorindicator','residual','type','powersubspace','restart',3);
[qrad,result]=solve(GSD,Ksto,fsto);

%% 
Kstoc = (Ksto);
fstoc = mat2cell(fsto);
GSD = GSDSOLVER('tol',1e-4,'nbfoncmax',6,...
    'display',true,'update',true,'errorindicator','residual','type','arnoldi','restart',3);
[qrad,result]=solve(GSD,Kstoc,fstoc);

GSD = GSDSOLVER('tol',1e-4,'nbfoncmax',6,...
    'display',true,'update',true,'errorindicator','residual','type','power','restart',3);
[qrad,result]=solve(GSD,Kstoc,fstoc);

GSD = GSDSOLVER('tol',1e-4,'nbfoncmax',2,...
    'display',true,'update',true,'errorindicator','residual','type','powersubspace','restart',3);
[qrad,result]=solve(GSD,Kstoc,fstoc);

