function [Xc,L] = PCidentification(PC,Xs,varargin)
% function [Xc,L] = PCidentification(PC,Xs)
% identification d'unn vecteur aleatoire sur le chaos polynomial
% PC : POLYCHAOS
% Xs : n-by-N array (N realisations d'un vecteur aleatoire de taille n)
%
% function Xc = PCidentification(PC,Xs,'nbsestim',nbsestim)
% 'nbsestim' : nombre de realisations utilisees pour l'estimation des pdf 
% marginales dans l'estimation de la vraisemblance approchee
%
% function Xc = PCidentification(PC,Xs,'parametrized')
% Variete de Stiefel parametree
%
% function Xc = PCidentification(PC,Xs,'nocenter')
% On ne determine pas la moyenne a priori (on rajoute donc n coefficients a identifier)
% (Permet de traiter les cas o� les polynomes du chaos ne sont pas tous � 
% esp�rance nulle)
%
% function Xc = PCidentification(PC,Xs,'initialguess',Xc0)
% Xc0 decomposition niitiale sur le chaos PC : initialisation de l'algo a
% gr
%
% ---------------------------------------
% CHOIX DE L'OPTIMISATION
% ----------------------------------------
%
% function Xc = PCidentification(PC,Xs,'rs',nbrandini)
% 'rs' : utilisation de random search
% nbrandini : nombre de tirage du random search optimization
% random search peut servir d'initialisation a un autre algo d'optimisation
%
% function Xc = PCidentification(PC,Xs,optimalgo,nbiter)
% optimalgo = condor ou fmin
% nbiter : nombre d'iterations 
%
% ---------------------------------------
% AFFICHAGE
% ----------------------------------------
%
% function Xc = PCidentification(PC,Xs,'illustrate')
% affichage de la fonction des iteres
%
% function Xc = PCidentification(PC,Xs,'plotfun',nbpoints)
% affichage de la fonction de vraisemblance (nbpoints : nombre de points par dimension)
% 
% ----------------------------------------
% EXPLICATION DE LA METHODE UTILISEE      
% Utilisation du maximum de vraisemblance
% ----------------------------------------
% Le vecteur aleatoire est transforme pour avoir une moyenne nulle et 
% une matrice de correlation correlation identite.
% La recherche des coefficients du chaos (n-by-(P+1) coefficients) 
% se ramene donc a la recherche 
% d'une matrice orthogonale n-by-P not�e A : on a donc un probleme d'optimisation 
% sur une variete de Stiefel compacte). 
%
%       max(L(A)) avec A matrice orthogonale
%
% Soit X le vecteur aleatore qu'on cherche a exprimer sur le chaos. A sont
% les coefficient du chaos. Soit Xs(:,k) les realisations du vecteur. La
% vraisemblance s'ecrit
%   L(A) = prod_k pX(Xs(:,k);A) (pX : densite de proba de X)
% et la vraisemblance approchee 
%   L(A) = prod_k prod_i pXi(Xs(i,k);A) (pXi : densite de proba marginale de Xi)
%
% %%%  Vision probleme de maximisation sous contrainte (par defaut)
% max(L(A)) sous la contrainte A*A'=I  (maximisation sous contrainte)
% utilisation d'un algorithme d'optimisation sous contrainte (condor ou fminsearch)
% random search utilise un tirage aleatoire d'une matrice orthogonale
%
% %%% Vision probleme de maximisation sans contrainte
% function Xc = PCidentification(PC,Xs,'parametrized')
% Si on a une parametrisation de la variete de Stiefel S(n,P) avec
%  (n*P-n*(n+1)/2) parametres
% Par exemple, pour n=1, S(1,P) est une hypershere parametree par P-1
% angles t1,t2,...tP-1
% A = [cos(t1),sin(t1)cos(t2),sin(t1)sin(t2)cos(t3),...,sin(t1)sin(t2)...sin(tP-1)]
%

if ischarin('initialguess',varargin)
    Xc0 = getcharin('initialguess',varargin);
    PC = getPC(Xc0);
end



global generalstiefel

parametrized = ischarin('parametrized',varargin);
illustrate = ischarin('illustrate',varargin);
illustraters = ischarin('illustraters',varargin);
minbounds = ischarin('minbounds',varargin);
nocenter = ischarin('nocenter',varargin);
mublock =  getcharin('mublock',varargin);
generalstiefel = getcharin('generalstiefel',varargin,0);

N = size(Xs,2);
n = size(Xs,1);
if ~isempty(mublock) 
    if n>1
        error('ca ne marche qu''avec un variable � identifier')
    end
    nbblock=length(mublock);
    P = length(PC)-nbblock;
    indpc = getindices(PC);
    orderpc = getorder(PC);
    polfe = getpoly(PC,1);
    if ~isa(polfe,'POLYFE') || getparam(polfe,'n')~=nbblock 
        error('le premier polynome du chaos doit etre un POLYFE compatible avec les blocks')
    end
    if orderpc(1)~=0
        error('l''ordre du chaos doit etre 0 suivant la dimension du POLYFE') 
    end
elseif ~nocenter
    P = length(PC)-1;
else
    P = length(PC);
    Hm = one(PC);
    if ~isempty(find(Hm(2:end)))
        error('mean value of chaos polynomials is not zero ! Utiliser argument ''nocenter''')
    end

end
dimvariete = n*P-n*(n+1)/2;


%% definition de fonctions et de la la fonction a minimiser

if ~nocenter & isempty(mublock)
    mu = sum(Xs,2)/N;
    Xtildes = (Xs - repmat(mu,1,N));
else
    Xtildes = Xs ;     
end

if ~ischarin('nosvd',varargin)
    [U,S,V] = svd(Xtildes) ;
    s = diag(S);
    rep = find(s/max(s)>1e-12);
    if length(rep)<n
        fprintf('\n !!! variables inter-d�pendantes ')
        fprintf('\n --> r�duction du nombre de variables al�atoires � identifier')    
        n = length(rep);
        U = U(:,rep);
        S = S(rep,rep);
        V=V(:,rep);

        [Xc,L] = PCidentification(PC,S*V',varargin{:}); 
        if ~nocenter & isempty(mublock)
            Xc = mu + U*Xc;
        else
            Xc = U*Xc;    
        end
        return
    end
end


fprintf('\n------------------------------------------------------')
fprintf('\n  IDENTIFICATION VARIABLES ALEATOIRES SUR LE CHAOS ')
fprintf('\n------------------------------------------------------')
fprintf('\nNombres de variables � identifier : %d',n)
fprintf('\nNombres d''echantillons            : %d',N)
fprintf('\ndimension du chaos                : %d',getM(PC))
fprintf('\ndegre du chaos                    : %d',max(getorder(PC)))
if parametrized
    fprintf('\nVision parametree de la variete de Stiefel')
    fprintf('\nnombre de parametres � identifier : %d',dimvariete)
else
    fprintf('\nVision non parametree de la variete de Stiefel -> optimisation sous contrainte')
    fprintf('\nnombre de parametres � identifier : %d',n*P)
    fprintf('\nnombre de contraintes             : %d',n*(n+1)/2)
end
fprintf('\n------------------------------------------------------\n')



if isempty(mublock)
    cov = (Xtildes*Xtildes')/(N-1);
    B = chol(cov);
    if ~nocenter
        afun = @(atilde) [mu,B'*atilde];
    else
        afun = @(atilde) [B'*atilde];    
    end

else
    cov = (Xs*Xs')/(N-1);
    mu = mean(Xs);
    indpc=getindices(PC);
    repm = find(indpc(:,2)==0);
    norepm = setdiff(1:length(PC),repm);
    mmm = getparam(polfe,'dx');
    firstcoeff = (mublock(:).*sqrt(mmm(:)))';
    B = cov-sum(firstcoeff.^2)
    if B<=eps
        error('mubllock conduit moment d''ordre 2 trop grand par rapport aux echantillons')   
    end
    B = sqrt(B);
    afun = @(atilde) goodplace(firstcoeff,B*atilde,repm,norepm);
end

atildefun = @(phi) eval_ortho_matrix_param(n,P,phi);
Xcfun = @(a) PCMATRIX(a,[n,1],PC);
Xcfunphi = @(phi) Xcfun(afun(atildefun(phi)));



Nbs = getcharin('nbsestim',varargin,1e4);
Nbsrs = getcharin('nbsestimrs',varargin,Nbs);

if dimvariete==0
    if n==1
        atilde = 1 ; 
    else
        error('indentification impossible : augmenter la dimension du chaos')    
    end
    a = afun(atilde);
    Xc = Xcfun(a);
    L = -likelihood(Xc,Xs,Nbs);
    compdfplot(Xs,{Xc},5e5,11,{'final'})
    return
end

if parametrized
    myfunvrais = @(phi,N) -likelihood(Xcfunphi(phi),Xs,N);
    [arg.lb,arg.ub] = ortho_matrix_param_bounds(n,P,varargin{:});
else
    myfunvrais = @(A,N) -likelihood(Xcfun(afun(reshape(A,n,P))),Xs,N);    
    arg.lb = -ones(n*P,1);
    arg.ub = ones(n*P,1);
end


%% affichage de la vrais

if ischarin('plotfun',varargin)
    npts = getcharin('plotfun',varargin);
    if (parametrized & dimvariete<3) | (~parametrized & (n*P<=2))
        for i=1:1
            plot_fun(@(phi) myfunvrais(phi,Nbs),arg.lb,arg.ub,npts);
            hold on
        end    
        pause(.5)
    end
end



%% random search  pour l'initialisation
nbrandini = getcharin('rs',varargin,0);

if nbrandini>0
    fprintf('\n-----------------------')
    fprintf('\n    RANDOM SEARCH    ')
    fprintf('\n-----------------------\n')
    fprintf('nombre d''iterations = %d\n',nbrandini)


    if parametrized 
        if ~ischarin('lhsrandomsearch',varargin)
            [phi0,L]= my_random_search(@() rand_ortho_matrix_param(n,P,arg.lb,arg.ub),...
                @(x) myfunvrais(x,Nbsrs),nbrandini,illustraters);
        else     
            [phi0,L]= my_lhsrandom_search(arg.lb,arg.ub,...
                @(x) myfunvrais(x,Nbsrs),nbrandini,illustraters);
        end

        atilde0 = atildefun(phi0);
        L = myfunvrais(phi0,Nbs);
    else
        [atilde0,L]= my_random_search(@() rand_ortho_matrix(n,P),...
            @(x) myfunvrais(x,Nbsrs),nbrandini,illustraters);  
        L = myfunvrais(atilde0,Nbs);
    end

    fprintf('\nvraisemblance apres random search  = %f\n',L)
else
    fprintf('\n-----------------------')
    fprintf('\n INITIALISATION SIMPLE ')
    fprintf('\n-----------------------\n')


    if ischarin('initialguess',varargin)
        Xc0 = getcharin('initialguess',varargin);

        if ~nocenter 
            Xc0d = double(Xc0-mean(Xc0));
            atilde0 = (B')\Xc0d(:,2:end); 
        else
            atilde0 = (B')\double(Xc0);    
        end
        if parametrized
            phi0 = eval_ortho_matrix_param_inverse(n,P,atilde0);
            L = myfunvrais(phi0,Nbs);
        else
            L = myfunvrais(atilde0,Nbs);
        end

    else
        if parametrized 
            phi0=zeros(dimvariete,1); 
            L= myfunvrais(phi0,Nbs);
            atilde0 = atildefun(phi0);
        else
            atilde0 = zeros(n,P);for i=n;atilde0(i,i)=1;end
            L = myfunvrais(atilde0,Nbs);
        end
    end
    fprintf('\nvraisemblance apres initialisation = %f\n',L)

end
a0 = afun(atilde0);   
Xc0 = Xcfun(a0);
Xc=Xc0;



%%
%compdfplot(Xs,{Xc0},5e5,11,{'initial guess'})


%% ALGO CONDOR
nbiter = getcharin('condor',varargin,0);

if nbiter>0 
    fprintf('\n-----------------------')
    fprintf('\n  ALGORITHME CONDOR ')
    fprintf('\n-----------------------\n')
    fprintf('nombres d''iterations = %d\n',nbiter)

    if parametrized
        arg.xstart=phi0;
        arg.f='myfun_condor';
        opt.myfun = @(x) myfunvrais(x,Nbs);

        [phi,L,lambda,chemin] = matlabcondor(0.3,0.01,nbiter,arg,opt);
        phi=phi';
        Xc = Xcfunphi(phi);

    else
        arg.xstart=reshape(atilde0,n*P,1);
        arg.f='myfun_condor';
        arg.nNLConstraints=n*(n+1) ; 
        opt.myfun = @(x) myfunvrais(x,Nbs);
        opt.mycon = @(o,J,x) myconpourcondor(o,J,x,n,P);
        arg.c = 'mycon_condor';

        [atilde,L,lambda,chemin] = matlabcondor(0.3,0.01,nbiter,arg,opt);
        atilde = reshape(atilde,n,P);
        a = afun(atilde);
        Xc = Xcfun(a);
%fprintf('verification de la contrainte')
%atilde*atilde'
%error('programmer condor en nonlineaire')        
    end
    Xccondor=Xc;
    fprintf('\nvraisemblance apres    condor     = %f\n',L)

    if illustrate
        for i=1:size(chemin,1)   
            st.fval = chemin(i,end-1);
            outfun(chemin(i,1:end-2),st,'iter','rs');
            pause(0.2)
        end
    end    
    compdfplot(Xs,{Xc0,Xccondor},5e5,11,{'initial guess','final condor'})
end

%% Algo Optimisation matlab

nbiter = getcharin('fmin',varargin,0);
if nbiter>0    
    fprintf('\n-----------------------')
    fprintf('\n  ALGORITHME FMIN ')
    fprintf('\n-----------------------\n')
    fprintf('nombres d''iterations = %d\n',nbiter)

    options=optimset;
    options = optimset(options,'Display','iter');
    options = optimset(options,'TolFun',1e-2);
    if parametrized
        tolX = max(abs(arg.ub-arg.lb))*1e-2;
        options = optimset(options,'TolX',tolX);
    else
        tolX = 1e-2;
    end
    options = optimset(options,'MaxIter',nbiter);
    options = optimset(options,'LargeScale','off');
    objfun = @(x) myfunvrais(x,Nbs);
    if illustrate
        options = optimset(options,'OutputFcn', @(a,b,c) outfun(a,b,c,'ro',objfun));
    end


%options = optimset(options,'MaxFunEvals',nbiter);
    if parametrized 
        if minbounds    
            phi = fmincon(objfun,phi0,[],[],[],[],arg.lb,arg.ub,[],options);
        else
            phi = fminsearch(objfun,phi0,options);
        end

        L = myfunvrais(phi,Nbs);
        options = optimset(options,'GradConstr','off');
        atilde = atildefun(phi);
    else
        options = optimset(options,'GradConstr','on');

        if minbounds    
            atilde =  fmincon(objfun,atilde0(:),...
                [],[],[],[],arg.lb,arg.ub ,@(x) mycon(x,n,P),options);
        else
            atilde =  fmincon(objfun,atilde0(:),...
                [],[],[],[],[],[] ,@(x) mycon(x,n,P),options);    
        end
        atilde = reshape(atilde,n,P);
        L = objfun(atilde);
%atilde =  fminsearch(@(x) myfunvrais(x,Nbs),atilde0(:),options);
%atilde = reshape(atilde,n,P);
    end
%fprintf('verification de la contrainte')
%atilde*atilde'
    a=afun(atilde);
    Xc = Xcfun(a);
    Xcfmin=Xc;
    fprintf('\nvraisemblance apres   fmin     = %f\n',L)

    if illustrate
        compdfplot(Xs,{Xc0,Xcfmin},5e5,11,{'initial guess','final fmin'});
    end

end

fprintf('\n------------------------------------------------------')
fprintf('\n                FIN DE L''IDENTIFICATION          ')
fprintf('\n------------------------------------------------------\n')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FONCTIONS ADDITIONNELLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,L,L0]=my_random_search(argfun,minfun,n,illustrate)

L = Inf;

for i=1:n
    pourcentage(i,n)
    X0 = argfun();
    L0(i) = minfun(X0);
    s.fval = L0(i);
    marker = 'kx';
    if L0(i) < L
        X = X0;
        L = L0(i);
        marker ='gx';
    end
    state = 'iter';
    if illustrate
        outfun(X0,s,state,marker,minfun)  ;
    end
end

return 

function [X,L,L0]=my_lhsrandom_search(lb,ub,minfun,n,illustrate)

r = lhsdesign(n,length(lb));
lb = repmat(lb(:)',n,1);
ub = repmat(ub(:)',n,1);
r = lb+r.*(ub-lb);

L = Inf;

for i=1:n
    pourcentage(i,n)

    X0 = r(i,:);
    L0(i) = minfun(X0);
    s.fval = L0(i);
    marker = 'kx';
    if L0(i) < L
        X = X0;
        L = L0(i);
        marker ='gx';
    end
    state = 'iter';
    if illustrate
        outfun(X0,s,state,marker,minfun)  ;
    end
end


return

function plot_fun(fun,lb,ub,n)
fprintf('\nAffichage de la fonction\n')
figure(100)
if length(lb)==2
    [X1,X2] = meshgrid(linspace(lb(1),ub(1),n),linspace(lb(2),ub(2),n)); 
    [n1,n2]=size(X1);
    f=zeros(n1,n2);
    for i=1:n1,
        for j=1:n2,
            pourcentage((i-1)*n2+j,n1*n2,33)
            f(i,j)=fun([X1(i,j);X2(i,j)]);
        end
    end
    mxf=max(max(f));
    mnf=min(min(f));
    df=mnf+(mxf-mnf)*(2.^(([0:20]/20).^2)-1);
    clf
    [v,h]=contour(X1,X2,f,df);
%clabel(v,h);  
elseif length(ub)==1
    X1 = linspace(lb(1),ub(1),n) ;
    n1 = length(X1);
    f=zeros(n1,1);  
    for i=1:n1,
        pourcentage(i,n1)
        f(i)=fun(X1(i));
    end
    clf
    plot(X1,f);
else
    error('pas programme en dimension superieure a 2')
end
return


function stop=outfun(x, s, state,marker,fun)
if strcmp(state,'iter')
    figure(100)
    hold on
    if nargin<4 | isempty(marker)
        marker='ro'; 
    end

    if length(x)==1
        if nargin<5 
            L = s.fval;
        else
            L = fun(x);
        end

        plot(x,L,marker,'markersize',10)    
    elseif length(x)==2
        plot(x(1),x(2),marker,'markersize',10)        
    else
%error('pas d''affichage possible')
    end
    pause(.1)
end
stop=0;
return

function compdfplot(Xs,Xcs,nbs,fig,leg)
nbs = ceil((nbs).^(1/size(Xs,1)));
figure(fig)
clf
if size(Xs,1)>1
    subplot(1,length(Xcs)+1,1)
    multipdfsampleplot(Xs,'npts',40);
    ax0=axis;
    cax0=caxis;
    title('pdf sample')
else
    pdfsampleplot(Xs,getcourbestyles(1,'nomarker')); 
end
hold on
leg = [{'pdf sample'} leg];
for i=1:length(Xcs)
    if size(Xs,1)>1
        subplot(1,length(Xcs)+1,1+i)
        multipdfplot(Xcs{i},'nbs',nbs,'npts',30);
        axis(ax0);
        caxis(cax0);
        title(leg{i+1});
    else    
        pdfplot(Xcs{i},getcourbestyles(i+1,'nomarker'),'nbs',nbs)
    end
end

if size(Xs,1)==1
    legend(leg{:})
    title('-log(likelihood)')
end
pause(.5)
return



function [c,ceq,GC,GCeq] = mycon(x,n,P)

c=[];
x = reshape(x,n,P);
k = 0;
ceq = zeros(n*(n+1)/2,1);
for i=1:n
    for j=1:i
        k=k+1;        
        ceq(k)=x(i,:)*x(j,:)'-(i==j);  
    end
end

if nargout>2
    GC = zeros(n*P,0);
    GCeq = zeros(n,P,n*(n+1)/2);
    k = 0 ; 
    for i=1:n
        for j=1:i
            k=k+1;           
            GCeq(i,:,k)=x(j,:);
            GCeq(j,:,k)=GCeq(j,:,k)+x(i,:);
        end
    end
    GCeq = reshape(GCeq,n*P,n*(n+1)/2);
end

return


function [output,error] = myconpourcondor(grad,J,x,n,P)
error=0;

if J>n*(n+1)/2
    fact=-1;
    J=mod(J-1,n*(n+1)/2)+1;
else
    fact=1;
end

x = reshape(x,n,P);
if grad==0
    k = 0;
    for i=1:n
        for j=1:i
            k=k+1;  
            if k==J
                output=fact*(x(i,:)*x(j,:)'-(i==j));  
                break
            end
        end
    end

else
    output = zeros(n,P);
    k = 0 ; 
    for i=1:n
        for j=1:i
            k=k+1; 
            if k==J
                output(i,:)=fact*x(j,:);
                output(j,:)=output(j,:)+fact*x(i,:);
                break
            end
        end
    end
    output = reshape(output,n*P,1);
end

return

function a = goodplace(a1,a2,rep1,rep2)

a=zeros(1,length(rep1)+length(rep2));
a(rep1)=a1;
a(rep2)=a2;

return
