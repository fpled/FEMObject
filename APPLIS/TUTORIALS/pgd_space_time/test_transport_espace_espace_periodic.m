
%vx = @(x) (x(:,2)-1/2);
%vy = @(x) (1/2-x(:,1));

rx = 100;ry=rx;
Sx = mesh(DOMAIN(1),rx);
Sx = createddlnode(Sx,DDL('u'));
Sx = addclperiodic(Sx,POINT(0),POINT(1),'u');
Sy = mesh(DOMAIN(1),ry);Sy = createddlnode(Sy,DDL('u'));
Sy = addclperiodic(Sy,POINT(0),POINT(1),'u');
D=DOMAIN(2);
[S,X,Y]=mesh(D,rx,ry);

S = createddlnode(S,DDL('u'));
S = addclperiodic(S,getedge(D,1),getedge(D,3),'u');
S = addclperiodic(S,getedge(D,4),getedge(D,2),'u');

x = getcoord(getnode(Sx));
y = getcoord(getnode(Sy));
%cx = y(:)-1/2;
%cy = 1/2-x(:);
cx = 1;
cy = 2.5;

% cx = py(y), cy=px(x)
Ax = calc_matrix(BILINFORM(1,1),Sx);
Ay = calc_matrix(BILINFORM(1,1),Sy);
Mx = calc_matrix(BILINFORM(0,0),Sx);
My = calc_matrix(BILINFORM(0,0),Sy);
Mpx = calc_matrix(BILINFORM(0,0,cy,0),Sx);
Mpy = calc_matrix(BILINFORM(0,0,cx,0),Sy);
Cx = calc_matrix(BILINFORM(0,1),Sx);
Cy = calc_matrix(BILINFORM(0,1),Sy);

nu = 1;
si = 10;
c = 25;
Asep = SEPMATRIX({Ax,My;Mx,Ay;Cx,Mpy;Mpx,Cy;Mx,My},[nu,nu,c,c,si]);

%
fx = LINFORM(0,exp(-(x(:)-1/4).^2/0.1^2*3),0);
fx = fx{Sx}(:);
fy = LINFORM(0,exp(-(y(:)-1/2).^2/0.1^2*3),0);
fy = fy{Sy}(:);
bsep = SEPMATRIX({fx,fy});
%

bb = expand(bsep);
bb = unfreevector(Sy,unfreevector(Sx,bb,2),1)
figure(3)
clf
surf(X,Y,bb'+1e-5,'EdgeColor','none')
axis off
axis square
view(2)


%%
uref=[];
Asep = setalpha(Asep,[nu,nu,50,50,si]);
adjoint =0;
resi = 1;
update = 0;
mmax=100;
alphaupdate=0;
%uref = [];
PGD = SEPSOLVER(getdim(Asep),'adjoint',adjoint,'maxorder',mmax,...
    'residual',resi,'maxiter',7,'update',update,'updateadjoint',false,...
    'itercrit',1e-2,'alphaupdate',alphaupdate,'errorindicator','residual',...
    'updatedim',1:getdim(Asep),'ortho',true,'inittype','one',...
    'reference',uref,...
    'righthandSD',1)%,'metric',Ametric);

[usep,res] = solve(Asep,bsep,PGD);

if adjoint
    col = 'ms';
    if update 
        resadu = res;
        usepadu = usep;
    else
        resad = res;
        usepad = usep;    
    end
elseif resi
    col = 'r^';
    if update
        resresiu = res;
        usepresiu = usep;
    else
        resresi = res;
        usepresi = usep;    
    end
else
    col = 'bo';
    if update
        resgalu = res;
        usepgalu = usep;
    else
        resgal = res;
        usepgal = usep; 
    end
end
if update
    col=[col '--'];
else
    col = [col '-'];
end

figure(23)
semilogy(res.error,col)
hold on

%
uplot = usep;
bb = expand(uplot);
bb = unfreevector(Sy,unfreevector(Sx,bb,2),1);

figure(3)
clf
surf(X,Y,bb'+1e-5,'EdgeColor','none')
axis off
axis square
view(2)

