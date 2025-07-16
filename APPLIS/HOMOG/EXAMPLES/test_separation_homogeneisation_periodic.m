
rx = 100;ry=rx;
Sx = mesh(DOMAIN(1),rx);
Sx = createddlnode(Sx,DDL('u'));
Sx = addclperiodic(Sx,POINT(0),POINT(1),'u');
Sy = mesh(DOMAIN(1),ry);
Sy = createddlnode(Sy,DDL('u'));
Sy = addclperiodic(Sy,POINT(0),POINT(1),'u');
D=DOMAIN(2);
[S,X,Y]=mesh(D,rx,ry);

S = createddlnode(S,DDL('u'));
S = addclperiodic(S,getedge(D,1),getedge(D,3),'u');
S = addclperiodic(S,getedge(D,4),getedge(D,2),'u');

x = getcoord(getnode(Sx));
y = getcoord(getnode(Sy));

% cx = py(y), cy=px(x)
a = setfree(BILINFORM(1,1),0);
m = setfree(BILINFORM(0,0),0);
nu = 1;
Ax = calc_matrix(a,Sx);
Ay = calc_matrix(a,Sy);
Mx = calc_matrix(m,Sx);
My = calc_matrix(m,Sy);
Asepfree = SEPMATRIX({Ax,My;Mx,Ay},[nu,nu]);

%%% ajout d'une heterogeneite %%
if 0
    ax = setfree(BILINFORM(1,1,exp(-20*(x-1.3*mean(x)).^2),0),0);
    ay = setfree(BILINFORM(1,1,exp(-20*(y-mean(y)).^2),0),0);
    mx = setfree(BILINFORM(0,0,exp(-20*(x-1.3*mean(x)).^2),0),0);
    my = setfree(BILINFORM(0,0,exp(-20*(y-mean(y)).^2),0),0);
    Axx = calc_matrix(ax,Sx);
    Ayy = calc_matrix(ay,Sy);
    Mxx = calc_matrix(mx,Sx);
    Myy = calc_matrix(my,Sy);
    Asepfree = Asepfree + SEPMATRIX({Axx,Myy;Mxx,Ayy},[nu,nu]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%  ajout d'une penalisation (inutile) de la moyenne
lx = setfree( LINFORM(0,1),0);
Lx = calc_vector(lx,Sx);
Ly = calc_vector(lx,Sy);
if 0
    epsi = 1e5;
    Asepfree = Asepfree + SEPMATRIX({Lx*Lx',Ly*Ly'},epsi);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
Asep = freematrix(Asepfree,1,Sy);
Asep = freematrix(Asep,2,Sx);

umacro = SEPMATRIX({x(:),ones(length(y),1)});
bmacro = -Asepfree*umacro;
bmacro = freevector(bmacro,2,Sx);
bmacro = freevector(bmacro,1,Sy);
bsep = bmacro ;
%% affichage du seconde membre
bb = expand(bsep);
bb = unfreevector(Sy,unfreevector(Sx,bb,2),1);
figure(3)
clf
surf(X,Y,bb','EdgeColor','none')
axis off
axis square
view(2)


%%


update = 0;
mmax=100;
errorindicator = 'none';
PGD = SEPSOLVER(getdim(Asep),'tol',1e-15,'maxorder',mmax,...
    'maxiter',7,'update',update,'updateadjoint',false,...
    'errorindicator',errorindicator,...
    'updatedim',1:getdim(Asep));

[usep,res] = solve(Asep,bsep,PGD);

%%%  moyenne et projection
u1 = expand(usep);
u1 = unfreevector(Sy,unfreevector(Sx,u1,2),1);
moy = Lx'*u1*Ly; % calcul de la moyenne
u1 = u1-moy;  % on projete sur l'espace des fonctions a moyenne nulle
fprintf('moyenne = %d \n',moy)

col = 'bo';
figure(23)
semilogy(res.error,col)
hold on
%

figure(3)
clf
surf(X,Y,u1,'EdgeColor','none')
axis off
axis square
view(2)


%%% un calcul particulier pris comme reference
if 1
    u2=u1;
else
    err=norm(u1-u2)/norm(u2);
    fprintf('erreur entre grandeurs projetees = %d \n',err)
end



