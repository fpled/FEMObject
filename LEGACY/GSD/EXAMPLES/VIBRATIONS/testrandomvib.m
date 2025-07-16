%% modele elements finis (S, f, K1, M1, ddloutput)
choixmesh = 1;

switch choixmesh
    case 1
        example_2Dsimple_gmsh
    case 2
        example_piece_gmsh
    case 3
        example_plaque_gmsh
    case 4
        example_trapeze_3D
    case 5
        M1 = 1;
        K1 = 1;
        f = 1;
        
        wband = [0,2];
        ddloutput = 1;
        alpha = @(w) 0.01;
        beta = @(w) 0.01;
end

%% DEFINITION DU MODELE ALEATOIRE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
isrand = 1;
nbvar = 1;
amort = 1;

if isrand
    %R = RANDVARS(RVLOGUNIFORM(0.1,1),nbvar);
    %R = RANDVARS(RVUNIFORM(.1,1),nbvar);
    R = RANDVARS(RVUNIFORM(0.5,1),nbvar);
    PCM = PCMODEL(R,'order',4,'fedim',1:nbvar,'femesh',repmat({10},1,nbvar));
    %PCM = PCMODEL(R,'order',8,'pcg')
    PC = getPC(PCM);
    E = PCM{1};
    %E = one(PC) ;
    %gamma = PCM{1};
    gamma = 1;
    if nbvar==1
        RHO = one(PC);
    else
        RHO = PCM{2};
    end
    
else
    E = 1;
    RHO = 1;
    gamma = 1;
    
end
M = RHO*M1;
K = E*K1;

if isrand
    if amort==0
        b = one(PC)*f;
        A = @(w) K - w^2*M;
        C = zero(PC)*sparse(size(K,1),size(M,1));
    else
        b = one(PC)*[f;0*f];
        Z = sparse(size(M1,1),size(M1,2));
        A = @(w) E* [K1,Z;Z,K1] + ...
            RHO*[-w^2*M1,Z;Z,-w^2*M1] + ...
            gamma*E*[Z,-w*beta(w)*K1;w*beta(w)*K1,Z] + ...
            gamma*RHO*[Z,-w*alpha(w)*M1;w*alpha(w)*M1,Z];
        Ar = @(w) E*K1-RHO*w^2*M1 ;
        Ai = @(w) (gamma*E)*w*beta(w)*K1 + (gamma*RHO)*w*alpha(w)*M1;
        C = @(w) gamma*(alpha(w)*M+beta(w)*K);
    end
elseif ~isrand
    if amort==0
        A = @(w) K-w^2*M ;
        b = f;
    else
        C = @(w) alpha(w)*M+beta(w)*K;
        A = @(w) [K-w^2*M,-gamma*w*C(w);gamma*w*C(w),K-w^2*M] ;
        b = [f;0*f];
    end
end

%% RESOLUTION CLASSIQUE
xi = random(RANDVARS(PC));
c = 1;
testr = 0 % pour tester une realisation
scanw = linspace(wband(1),wband(2),50);
%scanw = 0.23;
if isrand && testr==0
    urw = zeros(length(ddloutput),length(scanw),PC);
    uiw = zeros(length(ddloutput),length(scanw),PC);
else
    urw = zeros(length(ddloutput),length(scanw));
    uiw = zeros(length(ddloutput),length(scanw));
end

for i=1:length(scanw)
    %pourcentage(i,length(scanw))
    fprintf('%d / %d \n',i,length(scanw))
    w = scanw(i);
    
    if isrand
        if testr==1
            u = solve(randomeval(A(w),xi,RANDVARS(PC)),randomeval(b,xi,RANDVARS(PC)));
        else
            u = cgs(A(w),b,1e-5);
        end
    else
        u = solve(A(w),b);
    end
    ur = u(1:size(K,1));
    urw(:,i) = ur(ddloutput);
    if amort
        ui = u(size(K,1)+1:end);
        uiw(:,i) = ui(ddloutput);
    end
    
end

if ~israndom(urw)
    figure(100)
    hold on
    umod = sqrt(urw.^2+uiw.^2);
    plot(scanw,log(umod(c,:)));
end

urwref = urw;
uiwref = uiw;

%%
c = 1;
xi = random(RANDVARS(PC));
Kr = randomeval(K,xi);
Mr = randomeval(M,xi);
Cr = @(w) randomeval(C(w),xi);
[urw,uiw,umod,uphi] = calc_frf(scanw,Kr,Cr,Mr,f,ddloutput);
figure(101)
hold on
plot(scanw,log(umod(c,:)),'r');

%% RESOLUTION GSD
GSD = GSDSOLVER('tol',1e-5,'nbfoncmax',30,...
    'display',true,'update',true,'errorindicator','residual','type',...
    'power','restart',3,'orthocrit',1e-14,'nbfoncmaxsimul',1,...
    'inittype','random','pfixmax',3,'pfixtol',1e-3,'direct',true);
GSDA = setparam(GSD,'type','arnoldi');
GSDP = setparam(GSD,'type','power');

GSD = GSDP;
scanw = linspace(wband(1),wband(2),50);

urw = zeros(length(ddloutput),length(scanw),PC);
uiw = zeros(length(ddloutput),length(scanw),PC);
normresref = 0;
reuse = 1;

tol = 1e-2;
tolsolve = tol;
GSD = setparam(GSD,'tol',tolsolve);
GSD = setparam(GSD,'errorindicator','residual');

W = zeros(size(A(0),1),0);
m = size(W,2);
ajoutfun = zeros(1,length(scanw));
nbfun = zeros(1,length(scanw));
decompSD = 0;
decompSDtol = 1e-6;
for w=scanw
    fprintf('\n ----------------  w = %.3f  (%d/%d) ------------ \n',w,find(w==scanw),length(scanw))
    Aw = A(w);
    Aw = calc_ximasse(Aw,PC);
    % W = zeros(size(A(0),1),0);
    if size(W,2)>0 && reuse
        W = orth(full(W));
        lam = solve(W'*Aw*W,W'*b,tol/100);
        urad = PCRADIALMATRIX(W,[size(Aw,1),1],lam);
        
        if getm(urad)>12
            resi = spectral_decomposition(b-Aw*urad,'tol',tol/10);
        else
            resi = b-Aw*urad;
        end
        
        normresref = max(normresref,norm(b));
        
        epsi = norm(resi)/normresref;% norm(urad-u)/norm(u);
        fprintf('\n m = %d ; eps = %.2d\n',getm(urad),epsi)
    else
        epsi = 1;
    end
    
    if epsi>=tol
        if size(W,2)>0 && reuse
            GSD = setparam(GSD,'reuse',true);
            GSD = setparam(GSD,'update',true);
            [urad,result] = solve(GSD,Aw,b,urad);
        else
            GSD = setparam(GSD,'reuse',false);
            GSD = setparam(GSD,'update',true);
            [urad,result] = solve(GSD,Aw,b);
        end
        W = double(getV(urad));
        ajoutfun(w==scanw) = size(W,2)-m;
        m = size(W,2);
        %lam = solve(W'*A*W,W'*b);
        %urad = PCRADIALMATRIX(W,size(u),lam);
        %epsi = norm(b-A*urad)/norm(b);
        try
            fprintf('\n Ajout\n m = %d ; eps = %.2d\n',getm(urad),result.error(end))
        catch
            fprintf('\n Ajout\n m = %d ; eps < %.2d\n',getm(urad),tolsolve)
        end
        %norm(urad-u)/norm(u);
        if decompSD
            urad = spectral_decomposition(urad,'tol',decompSDtol)
            W = double(getV(urad));
            m = size(W,2);
        end
        
    end
    i = find(w==scanw);
    ur = urad(1:size(K,1));
    ui = urad(size(K,1)+1:2*size(K,1));
    urw(:,i) = expand(ur(ddloutput));
    uiw(:,i) = expand(ui(ddloutput));
    
    nbfun(w==scanw) = getm(urad);
    %pause(2)
end

urwrad = urw ;
uiwrad = uiw;

if ~reuse
    nbfunwithoutreuse = nbfun;
else
    nbfunwithreuse = nbfun;
    ajoutfunwithreuse = ajoutfun;
end

%%
figure(200)
clf
bar(scanw,[ajoutfunwithreuse;nbfunwithreuse-ajoutfunwithreuse]','stacked')
colormap(gray)
legend('Reuse - added functions','Reuse - Total number M')
%set(gca,'Ytick',[0:5:])
ddw = scanw(2)-scanw(1);
xlim([-ddw,max(scanw)+ddw])
% bar(scanw,ajoutfun,1)

hold on

%%
figure(200)
clf
if ~reuse
    bar(scanw,[nbfunwithoutreuse]','stacked','w')
else
    bar(scanw,[nbfunwithreuse]','stacked','w')
end
colormap(gray)
legend('Without Reuse')
%set(gca,'Ytick',[0:5:])
ddw = scanw(2)-scanw(1);
xlim([-ddw,max(scanw)+ddw])
% bar(scanw,ajoutfun,1)


%%
figure(1);
clf
%plot(Sb,'EdgeAlpha',1,'EdgeColor','none','FaceColor','b','FaceAlpha',1)
scanmode = [1,3:3:24];
Vr = unfreevector(S,W(1:size(K,1),:));
Vi = unfreevector(S,W(size(K,1)+1:end,:));
Vx = Vr(findddl(S,'UX'),:);
Vy = Vr(findddl(S,'UY'),:);
nbm = length(scanmode);
for k=scanmode
    fullsubplot(ceil(sqrt(nbm)),ceil(nbm/ceil(sqrt(nbm))),find(k==scanmode))
    text(1,-0.2,['| U_{' num2str(k) '}|' ],'FontSize',14)
    scal = 1/max(max(abs(V(:,k))))*max(max(abs(getcoord(getnode(S)))))*0.05;
    plot(FENODEFIELD(sqrt(Vy(:,k).^2+Vx(:,k).^2)),S+scal*V(:,k));
    cx=caxis();
    caxis([0,cx(2)*.8])
    pause(0.5)
    view(3)
end

%%
urw = urwref;
uiw = uiwref;
%%
urw = urwrad;
uiw = uiwrad;
%% AFFICHAGE DES FRF
c = 1
nbsimul = 1e5;
col='r';
if amort
    [urwr,xi] = random(urw(c,:),nbsimul,1);
    uiwr = randomeval(uiw(c,:),xi);
    uwmod = sqrt(urwr.^2+uiwr.^2);
else
    [urwr,xi] = random(urw(c,:),nbsimul,1);
    uwmod = abs(urwr);
end
Y = quantile(double(log(uwmod)),[0.001,.999],2)
%
figure(50)
hold on

monte = 1;
for kk=1:4
    if monte
        [Ks,xi] = random(K);
        Ms = randomeval(M,xi,RANDVARS(K));
        gammas = randomeval(gamma,xi);
        Cs = @(w) gammas*(alpha(w)*Ms+beta(w)*Ks);
        [urws,uiws,umods,uphis] = calc_frf(scanw,Ks,Cs,Ms,f,ddloutput(c));
    else
        [urr,xi] = random(urw(c,:));
        uir = randomeval(uiw(c,:),xi);
        umods = sqrt(urr.^2+uir.^2);
    end
    plot(scanw,log(umods))
    hold on
end
plotenveloppe(scanw,(Y),col,...
    'LineWidth',1,'FaceAlpha',0.5,'EdgeAlpha',1,'EdgeColor','k')

%%
figure(50)
ddw = scanw(2)-scanw(1);
xlim([-ddw,max(scanw)+ddw])

ax1 = gca;
ax2 = axes('Position',get(ax1,'Position'),...
    'XAxisLocation','bottom',...
    'YAxisLocation','right',...
    'Color','none',...
    'XColor','k','YColor','k')
hold on
bar(scanw,ajoutfun,'stacked','FaceColor','k')
ylim([0,6])
set(gca,'Ytick',[0:6])
xlim([-ddw,max(scanw)+ddw])

%colormap(gray)

%set(gca,'Ytick',[0:5:])


%%
w = 2.5;
Aw = A(w);
Aiw = Ai(w);
Arw = Ar(w);

rr = 1:size(Aw,1)/2;
ri = size(Aw,1)/2+1:size(Aw,1);

lr0 = rand(PC);lr0 = normalize(lr0);
li0 = rand(PC);li0 = normalize(li0);
l0 = [lr0;li0];
lr = lr0;
li = li0;
for kk=1:100
    Br = expect(Arw,lr,lr)+expect(Arw,li,li);
    Bi = expect(Aiw,lr,lr)+expect(Aiw,li,li);
    cr = expect(lr)*f;
    ci = expect(li)*f;
    U = [Br,-Bi;Bi,Br]\[cr;ci];
    U = U/norm(U);
    Ur = U(rr);
    Ui = U(ri);
    Br = Ur'*Arw*Ur+Ui'*Arw*Ui;
    Bi = Ur'*Aiw*Ur+Ui'*Aiw*Ui;
    cr = Ur'*f;
    ci = Ui'*f;
    l = solve([Br,-Bi;Bi,Br],one(PC)*[cr;ci]);
    lr = l(1);
    li = l(2);
    
    Ulr = Ur*lr-Ui*li ;
    Uli = Ui*lr+Ur*li;
    
    if kk>1
        DUlr = Ulr - Ulr0 ;
        DUli = Uli - Uli0 ;
        sqrt((norm(DUlr)^2+norm(DUli)^2)/(norm(Ulr)^2+norm(Uli)^2))
    end
    Ulr0 = Ulr;
    Uli0 = Uli;
end


%%
w = 2.5;
Aw = A(w);
r1 = 1:size(Aw,1)/2;
r2 = size(Aw,1)/2+1:size(Aw,1);

A11 = Aw(r1,r1);
A12 = Aw(r1,r2);
A21 = Aw(r2,r1);
A22 = Aw(r2,r2);
b1 = b(r1);
b2 = b(r2);

urad1 = PCRADIALMATRIX([length(r1),1],PC);
urad2 = PCRADIALMATRIX([length(r2),1],PC);
W1 = zeros(length(r1),0);
W2 = zeros(length(r2),0);
stoerror = [];
for mm=1:10
    l10 = rand(PC);l10 = normalize(l10);
    l20 = rand(PC);l20 = normalize(l20);
    l1 = l10;
    l2 = l20;
    for kk=1:2
        B11 = expect(A11,l1,l1);
        B12 = expect(A12,l1,l2);
        B21 = expect(A21,l2,l1);
        B22 = expect(A22,l2,l2);
        c1 = expect(b1,l1);
        c2 = expect(b2,l2);
        U = [B11,B12;B21,B22]\[c1;c2];
        U=U/norm(U);
        U1 = U(r1);
        U2 = U(r2);
        U1 = U1/norm(U1);
        U2 = U2/norm(U2);
        
        B11 = U1'*A11*U1;
        B12 = U1'*A12*U2;
        B21 = U2'*A21*U1;
        B22 = U2'*A22*U2;
        c1 = U1'*b1;
        c2 = U2'*b2;
        l = solve([B11,B12;B21,B22],[c1;c2]);
        l1 = l(1);
        l2 = l(2);
        Ul1 = U1*l1;
        Ul2 = U2*l2;
        if kk>1
            DUl1 = Ul1 - Ul10 ;
            DUl2 = Ul2 - Ul20 ;
            %sqrt((norm(DUl1)^2+norm(DUl2)^2)/(norm(Ul1)^2+norm(Ul2)^2))
            %norm(l1-l10)
        end
        Ul10 = Ul1;
        Ul20 = Ul2;
        l10 = l1;
        l20 = l2;
    end
    fprintf('\nnombre de fonction %d',mm)
    W1 = [W1,U1];
    W2 = [W2,U2];
    
    B11 = W1'*A11*W1;
    B12 = W1'*A12*W2;
    B21 = W2'*A21*W1;
    B22 = W2'*A22*W2;
    c1 = W1'*b(r1);
    c2 = W2'*b(r2);
    l = solve([B11,B12;B21,B22],[c1;c2]);
    l1 = l(1:mm);
    l2 = l(mm+1:2*mm);
    urad1 = PCRADIALMATRIX(W1,size(urad1),l1);
    urad2 = PCRADIALMATRIX(W2,size(urad2),l2);
    b1 = b(r1) - A11*urad1 - A12*urad2;
    b2 = b(r2) - A21*urad1 - A22*urad2;
    
    bres = [b1;b2];
    stoerror = [stoerror,norm(bres)/norm(b)];
    fprintf(' -- residu = %d',stoerror(end));
end

GSDP = setparam(GSDP,'nbfoncmax',mm);
[urad,result] = solve(GSDP,Aw,b);

