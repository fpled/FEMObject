
rx=1000;
rt=100;
dims = 1;
calcref = 1;
k=1e-1;

T = TIMEMODEL(0,1,rt);
N = DGTIMESOLVER(T,0);
Mt = getMmatrix(N);
Dt = getDmatrix(N,'basic');
nt = length(N);
ft = getMmatrix(N)*full(double(one(N))');

Sx = mesh(DOMAIN(1),rx);
Sx =  createddlnode(Sx,DDL('u'),DDL('g'));
Sx = addcl(Sx,POINT([0;1]),'u',0);
nx = getnbddlfree(Sx);

ax = BILINFORM(1,1);
cx = BILINFORM(0,1);
mx = BILINFORM(0,0);
Ax = ax{Sx}(:,:);
Mx = mx{Sx}(:,:);
Cx = cx{Sx}(:,:);
lx = LINFORM(0,1);
fx = lx{Sx}(:);

if dims==1
    S=Sx;
    M = Mx;
    A = k*Ax;
    fS=fx;
    MS=M;
    AS=A;
elseif dims==2
    Sy=Sx;
    ny = getnbddlfree(Sy);
    fy=fx;
    Ay = Ax;
    Cy = Cx;
    My=Mx;
    [S,X,Y] = mesh(DOMAIN(2),rx,rx);
    S =  createddlnode(S,DDL('u'),DDL('g'));
    S = addcl(S,[],'u',0);
    a = BILINFORM(1,1);
    c = BILINFORM(0,1,[1;1]);
    m = BILINFORM(0,0);
    ca=5;
    AS = k*a{S}(:,:)+ca*c{S}(:,:);
    MS = m{S}(:,:);
    l = LINFORM(0,1);
    fS = l{S}(:);

elseif dims==3
    Sy=Sx;
    Sz=Sx;
    ny = getnbddlfree(Sy);
    nz = getnbddlfree(Sz);
    fy=fx;
    fz=fx;

    Ay = Ax;
    My=Mx;
    Az = Ax;
    Mz=Mx;
    [S2,X2,Y2] = mesh(DOMAIN(2),rx,rx,rx);
    S2 =  createddlnode(S2,DDL('u'),DDL('g'));
    if calcref
        [S,X,Y,Z] = mesh(DOMAIN(3),rx,rx,rx);
        S =  createddlnode(S,DDL('u'),DDL('g'));
        S = addcl(S,[],'u',0);
        a = BILINFORM(1,1);
        m = BILINFORM(0,0);
        AS = k*a{S}(:,:);
        MS = m{S}(:,:);
        l = LINFORM(0,1);
        fS = l{S}(:);
    end
end

% reference 
if calcref
    N = setparam(N,'display',true)
    uref = dsolve(N,fS*one(N),MS,AS);
    uref = unfreevector(S,uref);

%
    if dims==1
        s=[size(uref,1),length(N)];
        solref = reshape(double(uref),s);
        solref = freevector(Sx,solref,1);
    elseif dims==2
        s=[size(X,1),size(X,2),length(N)];    
        solref = reshape(double(uref),s);
        solref = freevector(Sx,solref,1);
        solref = freevector(Sy,solref,2);
    elseif dims==3
        s=[size(X,1),size(X,2),size(X,3),length(N)];    
        solref = reshape(double(uref),s);
        solref = freevector(Sx,solref,1);
        solref = freevector(Sy,solref,2);
        solref = freevector(Sz,solref,3);
    end
end


%% affichage reference 2D
gg=reshape(solref,size(solref,1)*size(solref,2),size(solref,3));
figure(45)
clf
evol(TIMEMATRIX(gg,N),S)
%% SEPARE
if calcref==0
    solref=[];
end
maxorder = 30;
maxiter=5;
update = 1;
alphaupdate = 0;
updatedim = 2%1:dims+1;
adjoint =0;
tol=1e-17
updatestep=1;
itercritupdate=1e-10;
if adjoint
    col = 'r';
else
    col = 'b';    
end
if update
    col=[col,'-^'];
else
    col=[col,'-*'];
end

errorindicator = 'residual';


if dims==1
    A = SEPMATRIX({Mx,Dt;k*Ax,Mt});
    b = SEPMATRIX({fx,ft});
    solver = SEPSOLVER(getdim(A),'tol',tol,'alphaupdate',alphaupdate,...
        'update',update,'updatedim',updatedim,...
        'maxorder',maxorder,'maxiter',maxiter,'reference',solref,...
        'adjoint',adjoint,'updatestep',updatestep,...
        'itercritupdate',itercritupdate,...
        'errorindicator',errorindicator,'updateadjoint',false,'updateeps',0,'ortho',1)

    [u,result] = solve(A,b,solver);

    u = unfreevector(u,1,Sx);
    alpha = getalpha(u);
    Vx = extractvectors(u,1);
    Vt = extractvectors(u,2);

elseif dims==2

    A = SEPMATRIX({Mx,My,Dt;k*Ax,My,Mt;k*Mx,Ay,Mt});
    A = A + SEPMATRIX({ca*Cx,My,Mt;Mx,ca*Cy,Mt});
    b = SEPMATRIX({fx,fy,ft});

    solver = SEPSOLVER(getdim(A),'tol',tol,'alphaupdate',alphaupdate,...
        'update',update,'updatedim',updatedim,...
        'maxorder',maxorder,'maxiter',maxiter,'reference',solref,...
        'adjoint',adjoint,'updatestep',updatestep,...
        'itercritupdate',itercritupdate,...
        'errorindicator',errorindicator,'updateadjoint',false);
    [u,result] = solve(A,b,solver);

    u = unfreevector(unfreevector(u,2,Sy),1,Sx);
    alpha = getalpha(u);
    Vx = extractvectors(u,1);
    Vy = extractvectors(u,2);
    Vt = extractvectors(u,3);

elseif dims==3

    A = SEPMATRIX({Mx,My,Mz,Dt;k*Ax,My,Mz,Mt;Mx,k*Ay,Mz,Mt;Mx,My,k*Az,Mt});
    b = SEPMATRIX({fx,fy,fz,ft});

    solver = SEPSOLVER(getdim(A),'tol',tol,'alphaupdate',alphaupdate,...
        'update',update,'updatedim',updatedim,...
        'maxorder',maxorder,'maxiter',maxiter,'reference',solref,...
        'adjoint',adjoint,'updatestep',updatestep,...
        'itercritupdate',itercritupdate,...
        'errorindicator',errorindicator,'updateadjoint',false)
    [u,result] = solve(A,b,solver);

    alpha = getalpha(u);
    Vx = unfreevector(Sx,extractvectors(u,1));
    Vy = unfreevector(Sy,extractvectors(u,2));
    Vz = unfreevector(Sz,extractvectors(u,3));
    Vt = extractvectors(u,4);

end
figure(17);semilogy(result.error,col);
hold on

%%
col = 'y->'
[u,result] = solve(A'*A,A'*b,tol,'alphaupdate',alphaupdate,...
    'update',update,'updatedim',updatedim,...
    'maxorder',maxorder,'maxiter',maxiter,'reference',solref,...
    'adjoint',0);
figure(17);semilogy(result.error,col);
hold on


%% AFFICHE 1D

ugsd = Vx*diag(alpha)*Vt';
ugsd = TIMEMATRIX(ugsd,N,[size(ugsd,1),1]);

figure(5)
evol(ugsd,Sx,'courbe','rescale','y')

%% AFFICHE 2D
Vxy = zeros(numel(X),length(alpha));

for i=1:length(alpha)
    vtemp = Vx(:,i)*Vy(:,i)';    
    Vxy(:,i)=vtemp(:);
end

ugsd = Vxy*diag(alpha)*Vt';
ugsd = TIMEMATRIX(ugsd,N,[size(ugsd,1),1]);


figure(5)
evol(ugsd,S,'surface','rescale','z')

%% AFFICHE 3D coupe

coupe=13;
Vxyz = zeros(numel(X2),length(alpha));
for i=1:length(alpha)
    vtemp = Vx(:,i)*Vy(:,i)'*Vz(coupe,i);
    Vxyz(:,i)=vtemp(:);
end

ugsd2 = Vxyz*diag(alpha)*Vt';
ugsd2 = TIMEMATRIX(ugsd2,N,[size(ugsd2,1),1]);

figure(5)
evol(ugsd2,S2)%,'rescale','z','Color')



%% AFFICHE 3D

Vxyz = zeros(numel(X),length(alpha));    


for i=1:length(alpha)
    vtemp = Vx(:,i)*Vy(:,i)';
    vtemp = vtemp(:)*Vz(:,i)';    
    Vxyz(:,i)=vtemp(:);
end

ugsd = Vxyz*diag(alpha)*Vt';
ugsd = TIMEMATRIX(ugsd,N,[size(ugsd,1),1]);

figure(5)
evol(ugsd,keepelem(S,300:400),'colorbar','rescale','Color')    

