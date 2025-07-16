% exemple du PROCEEDING ECT 2010
%

%vx = @(x) (x(:,2)-1/2);
%vy = @(x) (1/2-x(:,1));

rx = 70;ry=rx;
Sx = mesh(DOMAIN(1),rx);
Sx = createddlnode(Sx,DDL('u'));
Sx = addclperiodic(Sx,POINT(0),POINT(1),'u');
Sy = mesh(DOMAIN(1),ry);Sy = createddlnode(Sy,DDL('u'));
Sy = addclperiodic(Sy,POINT(0),POINT(1),'u');


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


%S = createddlnode(S,DDL('u'));
%S = addclperiodic(S,getedge(D,2),getedge(D,4),'u');
%S = addclperiodic(S,getedge(D,1),getedge(D,3),'u');
D=DOMAIN(2);
[S,X,Y]=mesh(D,rx,ry);
mat = FOUR_ISOT('k',1,'b',[cx;cy],'c',1);
S = setmaterial(S,mat);
S=final(S);
S = addclperiodic(S,getedge(D,2),getedge(D,4),'T');
S = addclperiodic(S,getedge(D,1),getedge(D,3),'T');
M = calc_freematrix(S,@mass);
A = calc_freematrix(S,@diff);
mat = FOUR_ISOT('k',1,'b',[cx;0],'c',1);S = setmaterial(S,mat);
Ccx = calc_freematrix(S,@adv);
mat = FOUR_ISOT('k',1,'b',[0;cy],'c',1);S = setmaterial(S,mat);
Ccy = calc_freematrix(S,@adv);
fnode = @(x) exp(-(x(1)-1/4).^2/0.1^2*3).*exp(-(x(2)-1/2).^2/0.1^2*3);
f = bodyload(S,[],'QN',fnode);
%A =  calc_matrix(BILINFORM(1,1),S);
%C =  calc_matrix(BILINFORM(0,1,[cx;cy]),S);
%M = calc_matrix(BILINFORM(0,0),S);
%f = (LINFORM(0,exp(-(X(:)-1/4).^2/0.1^2*3).*exp(-(Y(:)-1/2).^2/0.1^2*3),0));
%f = f{S}(:);

nu = 1;
si = 10;
c = 50;

%
fx = LINFORM(0,exp(-(x(:)-1/4).^2/0.1^2*3),0);
fx = fx{Sx}(:);
fy = LINFORM(0,exp(-(y(:)-1/2).^2/0.1^2*3),0);
fy = fy{Sy}(:);

%
mmm=2;
if mmm
    RV = RANDVARS(RVLOGUNIFORM(1,10),RVUNIFORM(0,2.5));
    [RVpc,PC] = PCTPMODEL(RV,'fedim',1:length(RV),'femesh',repmat({14},1,length(RV)),'order',repmat(4,1,length(RV)));
    Mxione = get_ximasse(calc_ximasse(one(PC)));
    Mxinu = get_ximasse(RVpc{1});
    Mxicx = get_ximasse(RVpc{2});
else
    RV = RANDVARS(RVLOGUNIFORM(1,10));
    [RVpc,PC] = PCTPMODEL(RV,'fedim',1:length(RV),'femesh',repmat({14},1,length(RV)),'order',repmat(4,1,length(RV)));
    Mxione = get_ximasse(calc_ximasse(one(PC)));
    Mxinu = get_ximasse(RVpc{1});
    Mxicx = Mxione;
end
fxione = getphi(one(PC));
Asep = SEPMATRIX([{Ax,My},Mxinu;{Mx,Ay},Mxinu;{Cx,Mpy},Mxicx;{Mpx,Cy},Mxione;{Mx,My},Mxione]);
bsep = SEPMATRIX([{fx,fy},fxione]);
%Asep = SEPMATRIX({Ax,My,Mxinu;Mx,Ay,Mxinu;Cx,Mpy,Mxione;Mpx,Cy,Mxione;Mx,My,Mxione});
%bsep = SEPMATRIX({fx,fy,fxione});



%%

uref=[];
Asep = setalpha(Asep,[nu,nu,c,c,si]);

adjoint =0;
resi = 0;
update = 0;
if update
    mmax=50;
else
    mmax = 100;
end
alphaupdate=0;
%uref = [];
PGD = SEPSOLVER(getdim(Asep),'tol',1e-12,'adjoint',adjoint,'maxorder',mmax,...
    'residual',resi,'maxiter',5,'update',update,'updateadjoint',false,...
    'itercrit',1e-2,'alphaupdate',alphaupdate,'errorindicator','none',...
    'updatedim',1:getdim(Asep),'ortho',false,'inittype','one',...
    'reference',uref,'updatestep',1)%,'metric',Ametric);

[usep,res] = solve(Asep,bsep,PGD);

if adjoint
    col = 'gs';
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

figure(24)
semilogy(res.error,col)
hold on


%%
uplot = usep;
%H = polyval(PC,1/2);
[H,r] = random(PC);

factnu = randomeval(RVpc{1},r);
if mmm==2
    factcx = randomeval(RVpc{2},r);
    temp = mtimes(H{1},uplot,3);
    temp = mtimes(H{2},temp,4);
else
    temp = mtimes(H{1}',uplot,3);
end
usepr = expand(temp);
usepr = unfreevector(Sy,unfreevector(Sx,usepr,2),1)';
cax=caxis;
figure(3)
clf
%surf(X,Y,bb','EdgeColor','none')
plot(usepr(:),S)
axis off
axis square
view(2)

if mmm==2
    ur = (nu*factnu*A+c*(factcx*Ccx+Ccy)+si*M)\f;
else
    ur = (nu*factnu*A+c*(Ccx+Ccy)+si*M)\f;
end
ur = unfreevector(S,ur);
nn1=getnodenextto(getnode(S),POINT([0,1]));
nn2=getnodenextto(getnode(S),POINT([0,0]));
ur(getnumber(nn1))=ur(getnumber(nn2));

figure(4)
clf
%bb=reshape(ur,size(X,1),size(X,2));
%surf(X,Y,bb','EdgeColor','none')
plot(ur,S);
%caxis(cax);
axis off
axis square
view(2)

norm(usepr(:)-ur(:))/norm(ur(:))

%%
uplot = usep;
nn1=getnodenextto(getnode(S),POINT([0,1]));
nn2=getnodenextto(getnode(S),POINT([0,0]));
errnum = 0;
errden = 0;
for ii=1:30
    [H,r] = random(PC);
    factnu = randomeval(RVpc{1},r);
    temp = mtimes(H{1},uplot,3);
    if mmm==2
        temp = mtimes(H{2},temp,4);
        factcx = randomeval(RVpc{2},r);
    else
        factcx = 1;
    end

    usepr = expand(temp);
    usepr = unfreevector(Sy,unfreevector(Sx,usepr,2),1)';
    ur = (nu*factnu*A+c*(factcx*Ccx+Ccy)+si*M)\f;
    ur = unfreevector(S,ur);
    ur(getnumber(nn1))=ur(getnumber(nn2));
    errnum = errnum + norm(usepr(:)-ur(:))^2;
    errden = errden + norm(ur(:))^2;
end
err = sqrt(errnum/errden);
fprintf('\n L2 error = %d\n',err)


%% COMPARASION PGD et REFERENCE
for numsample = 6
    if mmm==2
        r = [0.5,0.5;0.1,0.1;0.1,0.9;0.9,0.1;0.9,0.9;1,1];
    else
        r = [0.1,0.3,0.5,0.7,0.9]';
    end
    r = r(numsample,:);
%[H,r] = random(PC);
    H = polyval(PC,r);
%
    factnu = randomeval(RVpc{1},r);
    temp = mtimes(H{1},uplot,3);
    if mmm==2
        temp = mtimes(H{2},temp,4);
        factcx = randomeval(RVpc{2},r);
    else
        factcx = 1;
    end

    if mmm==2
        ur = (nu*factnu*A+c*(factcx*Ccx+Ccy)+si*M)\f;
    else
        ur = (nu*factnu*A+c*(Ccx+Ccy)+si*M)\f;
    end

    ur = unfreevector(S,ur);
    ur(getnumber(nn1))=ur(getnumber(nn2));
    ur = reshape(ur,size(X,2),size(X,1));
%

    figure(33);
    clf;
    cax = [min(ur(:)),max(ur(:))];
    cval = linspace(cax(1)+(cax(2)-cax(1))*0.01,cax(2)-(cax(2)-cax(1))*0.01,20);
%surf(X,Y,ur,'EdgeColor','none');
%colormap('default')%colormap(flipud(bone))
    plot(DOMAIN(2),'LineWidth',2,'EdgeColor','k');
%plot(create_boundary(S))
%plot(DOMAIN(2),'FaceColor','w')
    hold on
    [Cref,Href]=contour(X,Y,ur,cval);
    view(2);axis image;axis off
    caxis(cax)

    myprint('C:/REDACTION/figs/ect2010_ex2/',['ex2_plot_ref_sample' num2str(numsample)],'epsc2')

%
    for cas=1:3
        for update=0:1
%cas = 1;update =1;
            if cas==3
                motsep =   'usepad';
                mot = 'minimax';
            elseif cas==2
                motsep =   'usepresi';
                mot = 'minres';
            elseif cas==1
                motsep = 'usepgal';
                mot = 'galerkin';
            end
            if update
                mot=['updated_' mot];
                motsep = [motsep 'u'];
                scanmm = [1:5,10:10:50];
            else
                mot=['progressive_' mot];
                scanmm = [1:5,10:20:60,100];
            end
            eval(['useptemp = ' motsep]);

            for mm = scanmm;
                temp = mtimes(H{1},truncate(useptemp,1:mm),3);
                if mmm==2
                    temp = mtimes(H{2},truncate(temp,1:mm),4);
                end
                temp=expand(temp);
                temp = unfreevector(Sy,unfreevector(Sx,temp,2),1)';

                figure(34);clf;
%surf(X(2:end-1,2:end-1),Y(2:end-1,2:end-1),ee,'EdgeColor','none');
                contour(X,Y,temp,cval);
                view(2);axis image;axis off
                plot(DOMAIN(2),'LineWidth',2);
                caxis(cax)
                view(2);axis image;axis off

                myprint('C:/REDACTION/figs/ect2010_ex2/',['ex2_plot_sample' num2str(numsample) '_' mot '_m' num2str(mm)],'epsc2')
            end
        end
    end

end

%% CONVERGENCE PGD progressive
%optionsplot = {'LineWidth',2,'MarkerSize',8};
optionsplot={};
imax = min(80,length(resgal.error));
errmin = 1e-6;
figure(24)
clf
%I = find(resref.error(1:imax)>errmin);
%semilogy(resref.error(1:I(end)),'k-*',optionsplot{:})
%hold on
I = find(resgal.error(1:imax)>errmin);
semilogy(I,resgal.error(I),'b-o',optionsplot{:})
hold on
I = find(resresi.error(1:imax)>errmin);
semilogy(I,resresi.error(I),'r-^',optionsplot{:})
I = find(resad.error(1:imax)>errmin);
semilogy(I,resad.error(I),'m-s',optionsplot{:})
I = find(resref.error(1:min(imax,end))>errmin);
semilogy(I,resref.error(I),'k-*',optionsplot{:})
%I = find(resrefuu.error(1:min(imax,end))>errmin);
%semilogy(I,resrefuu.error(I),'k--v',optionsplot{:})
ylim([1e-3,1])
set(gca,'FontSize',14)
legend('Galerkin-PGD','MinRes-PGD','Minimax-PGD','SVD','Location','northeast')

myprint('C:/REDACTION/figs/ect2010_ex2/','ex2_conv_progressive_PGD','epsc2')
%% CONVERGENCE PGD update
%optionsplot = {'LineWidth',2,'MarkerSize',8};
optionsplot={};
imax = min(40,length(resgal.error));
errmin = 1e-6;
figure(24)
clf
%I = find(resref.error(1:imax)>errmin);
%semilogy(resref.error(1:I(end)),'k-*',optionsplot{:})
%hold on
I = find(resgalu.error(1:imax)>errmin);
semilogy(I,resgalu.error(I),'b-o',optionsplot{:})
hold on
I = find(resresiu.error(1:imax)>errmin);
semilogy(I,resresiu.error(I),'r-^',optionsplot{:})
I = find(resadu.error(1:imax)>errmin);
semilogy(I,resadu.error(I),'m-s',optionsplot{:})
I = find(resrefu.error(1:imax)>errmin);
semilogy(I,resrefu.error(I),'k-*',optionsplot{:})
set(gca,'FontSize',14)
legend('Galerkin-PGD','MinRes-PGD','Minimax-PGD','SVD','Location','northeast')

myprint('C:/REDACTION/figs/ect2010_ex2/','ex2_conv_update_PGD','epsc2')

%% convergence avec multisvd

[usepref,resref]=multisvd(uref,'maxorder',100,'reference',uref);
[useprefu,resrefu]=multisvd(uref,'maxorder',50,'reference',uref,'update',1);
[useprefuu,resrefuu]=multisvd(uref,'maxorder',50,'reference',uref,'update',3);

figure(443)
semilogy(resref.error,'m-s',optionsplot{:})
hold on
semilogy(resrefu.error,'b-s',optionsplot{:},'MarkerFaceColor','b')
semilogy(resrefuu.error,'r-s',optionsplot{:},'MarkerFaceColor','r')


%% TEMPS DE CALCUL PROGRESSIVE
imax = min(90,length(resgal.error));
figure(44)
clf
rep = find(resgal.error(1:min(end,imax)));
semilogy(resgal.time(rep),resgal.error(rep),'b-o')
hold on
rep = find(resresi.error(1:min(end,imax)));
semilogy(resresi.time(rep),resresi.error(rep),'r-^')
rep = find(resad.error(1:min(end,imax)));
semilogy(resad.time(rep),resad.error(rep),'m-s')
xlim([0,200])
set(gca,'FontSize',14)
xlabel('(s)')
legend('Galerkin-PGD','MinRes-PGD','Minimax-PGD','Location','best')

myprint('C:/REDACTION/figs/ect2010_ex2/',['ex2_conv_time_progressive'],'epsc2')

%% TEMPS DE CALCUL UPDATED
imax = 240;
figure(44)
clf
rep = find(resgalu.error(1:min(end,imax)));
semilogy(resgalu.time(rep),resgalu.error(rep),'b-o')
hold on
rep = find(resresiu.error(1:min(end,imax)));
semilogy(resresiu.time(rep),resresiu.error(rep),'r-^')
rep = find(resadu.error(1:min(end,imax)));
semilogy(resadu.time(rep),resadu.error(rep),'m-s')
xlim([0,300])
set(gca,'FontSize',14)
xlabel('(s)')
legend('Galerkin-PGD','MinRes-PGD','Minimax-PGD','Location','northeast')

myprint('C:/REDACTION/figs/ect2010_ex2/',['ex2_conv_time_update'],'epsc2')

%% TEMPS DE CALCUL PROGRESSIVE ET UPDATE
figure(45)
clf
imax = 100;
rep = find(resgal.error(1:min(end,imax)));
semilogy(resgal.time(rep),resgal.error(rep),'b-o')
hold on
rep = find(resresi.error(1:min(end,imax)));
semilogy(resresi.time(rep),resresi.error(rep),'r-^')
rep = find(resad.error(1:min(end,imax)));
semilogy(resad.time(rep),resad.error(rep),'m-s')
%ylim([1e-10,1])
imax = 40;
rep = find(resgalu.error(1:min(end,imax)));
semilogy(resgalu.time(rep),resgalu.error(rep),'b-o','MarkerFaceColor','b')
hold on
rep = find(resresiu.error(1:min(end,imax)));
semilogy(resresiu.time(rep),resresiu.error(rep),'r-^','MarkerFaceColor','r')
rep = find(resadu.error(1:min(end,imax)));
semilogy(resadu.time(rep),resadu.error(rep),'m-s','MarkerFaceColor','m')
ylim([3e-3,1])
xlim([0,200])
set(gca,'FontSize',14)
xlabel('(s)')
legend('Galerkin-PGD','MinRes-PGD','Minimax-PGD',...
    'Galerkin-PGD-Updated','MinRes-PGD-Updated',...
    'Minimax-PGD-Updated','Location','southeast')

myprint('C:/REDACTION/figs/ect2010_ex2/',['ex2_conv_time_progressive_updated'],'epsc2')



%%

