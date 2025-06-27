function [u,result] = solve_power(GSD,A,b,ureuse,varargin)
% function u = solve_power(GSD,A,b,v,varargin)
% fonction solve GSD : resolution de Au=b
% A et b : PCMATRIX ou PCRADIALMATRIX
% v : pour la reutilisation (MULTIMATRIX ou PCRADIALMATRIX ou double)

dislay_ = getparam(GSD,'display');

fprintf('\n GSD solveur ... ')
if dislay_
    fprintf('\n')
end

param = getparam(GSD);
paramradial.nbfoncmax = getparam(GSD,'nbfoncmax');
paramradial.tol = getparam(GSD,'tol');
paramradial.pfixtol = getparam(GSD,'pfixtol');
paramradial.pfixmax = getparam(GSD,'pfixmax');
paramradial.reuse = getparam(GSD,'reuse');
paramradial.update = getparam(GSD,'update');
paramradial.updatefactor = getparam(GSD,'updatefactor');
paramradial.orthoduringpfix = getparam(GSD,'orthoduringpfix');
paramradial.adjoint =  getparam(GSD,'adjoint');

errorindicator = getparam(GSD,'errorindicator');
uref = getcharin('reference',varargin);

if ~isempty(uref)
    errorindicator='reference';
end

[A,b,PC]=pcsystemupdate(A,b,varargin{:});

if iscell(A) || iscell(b)
    % warning('on change le stockage des matrices aletaoires en double')
    % A = cell2mat(A);
    % b = cell2mat(b);
end

if getparam(GSD,'righthandSD');
    b = righthandSD(GSD,b);
end

P = getP(PC);
n=size(A,1);

errorpf=zeros(1,paramradial.pfixmax);
% result.error=zeros(1,paramradial.nbfoncmax);
% result.errorres=zeros(1,paramradial.nbfoncmax);
% result.errorrayg=zeros(1,paramradial.nbfoncmax);
% result.errorraygi=zeros(1,paramradial.nbfoncmax);
% result.time = zeros(1,paramradial.nbfoncmax);
% result.rayg=cell(1,paramradial.nbfoncmax);
% result.Rayg=cell(1,paramradial.nbfoncmax);
result.nfonctions = 0 ;
result.addedfunctions = 0;


l0 = zeros(1,PC);
l = zeros(1,PC);
U = sparse(n,1);
stoU = sparse(n,0);
U0 = sparse(n,1);
radalpha = [];


%% D�finition des normes

% Am = expect(A);
% Am = 1/2*(Am+Am');
Am = speye(n);

if ~strcmp(errorindicator,'rayleigh')
    mynorm=@(b) norm(b);
    normb = mynorm(b);
end

%% D�finition des solveurs
toliter = getparam(GSD,'toliter');
tolupdate = getparam(GSD,'tolupdate');

localstosolver = getparam(GSD,'localstosolver');
if ~isempty(localstosolver)
    localstosolveriter = @(A,b) solve(localstosolver,A,b,toliter);
    localstosolverupdate = @(A,b)  solve(localstosolver,A,b,tolupdate);
else
    directsolve = getparam(GSD,'direct');
    if directsolve
        localstosolveriter = @(aU,fU) solve(aU,fU);
        localstosolverupdate = @(aU,fU) solve(aU,fU);
    else
        localstosolveriter = @(aU,fU) cgs(aU,fU,toliter,[],'noupdate');
        localstosolverupdate = @(aU,fU) cgs(aU,fU,tolupdate,[],'noupdate');
    end
end

clock0 = clock;


%% R�utilisation
if (paramradial.reuse && nargin>=4 && ~isempty(ureuse))
    [u,flagreuse,resultini,bini]=solve_reuse(GSD,A,b,ureuse,varargin{:});
else
    flagreuse = true;
end
hasconv=false;

if flagreuse
    u = PCRADIALMATRIX([n,1],PC);
    bini = b;
    result.raygini=0;
    result.errorini=1;
else
    result.raygini=resultini.rayg;
    result.Raygini=resultini.Rayg;
    
    switch errorindicator
        case 'rayleigh'
            result.errorini=1;
        case 'reference'
            result.errorini=mynorm(u-uref)/mynorm(uref);
        case 'residual'
            result.errorini=mynorm(bini)/normb;
        case 'none'
            result.errorini=1;
    end
    if dislay_
        fprintf('  Initialisation : %3d fonctions',length(u))
    end
    if ~strcmp(errorindicator,'rayleigh')
        if dislay_
            fprintf(' - erreur  : %3d \n',result.errorini)
        end
        if result.errorini<getparam(GSD,'tolini')
            fprintf('  Convergence apres initialisation - erreur : %3d\n',result.errorini)
            hasconv=true;
            return
        end
    end
end

result.nfonctions = length(u);

bu = bini;

result.totaltime=etime(clock,clock0);

utilde = PCRADIALMATRIX([n,1],PC);

%% Construction des couples de fonctions
nbfoncmax = paramradial.nbfoncmax;

for j=1:nbfoncmax %%%%%%%%%%%%%%%% BOUCLE NOUVELLES FONCTIONS
    
    result.nfonctions = result.nfonctions + 1  ;
    
    switch getparam(GSD,'inittype')
        case 'one'
            l0 = one(PC);
        case 'random'
            l0 = rand(PC);
        case 'allone'
            l0 = ones(length(PC),1);
            l0 = PCMATRIX(l0,[1,1],PC);
            
    end
    l0tilde = l0;
    alpha0=norm(l0);l0=l0/alpha0;
    
    for i=1:paramradial.pfixmax %%%%%%%%%%%%%%%% BOUCLE POINT FIXE
        
        f0=expect(bu,l0tilde);
        L0 = expect(A,l0,l0tilde);
        
        perm = symrcm(L0);
        U(perm) = L0(perm,perm)\f0(perm);
        %rayg(i)=U'*f0;
        if paramradial.orthoduringpfix
            V = double(getV(u));
            if size(V,2)>=1
                U = U - V*(V'*Am*U);
            end
        end
        alpha = sqrt(U'*Am*U);
        U = U /alpha ;
        
        if param.saveiter
            result.lambda{j}{i}=l0;
        end
        
        if paramradial.adjoint
            f0=expect(l0,l0)*Am*U;
            L0=expect(A,l0,l0tilde)';
            Utilde = L0\f0;
            Utilde = Utilde/sqrt(Utilde'*Am*Utilde);
        else
            Utilde=U;
        end
        
        fU = Utilde'*bu;
        aU = Utilde'*A*U;
        
        [l,flag]=localstosolveriter(aU,fU);
        l = expand(l);
        
        
        if paramradial.adjoint
            
            fUtilde = (U'*Am*U)*l;
            aUtilde = Utilde'*A*U;
            
            [ltilde,flag]=localstosolveriter(aUtilde,fUtilde);
            ltilde = expand(ltilde);
            ltilde = ltilde/norm(ltilde);
        else
            ltilde = l;
        end
        
        result.rayg{j}(i) = full(expectmtimes(ltilde,cell2mat(fU)'));
        
        if param.saveiter
            result.U{j}{i}=U;
        end
        
        % errorpf(i)=double(sqrt(...
        %   (alpha0^2+alpha^2-2*alpha0*alpha*(U0'*Am*U)*(double(l0)*double(l)'))/...
        %                 (alpha^2)));
        errorpf(i)=abs((result.rayg{j}(i)-(i>1)*result.rayg{j}(i-(i>1)))/result.rayg{j}(i));
        
        % errorpf(i)=abs(alpha-alpha0)/(alpha+alpha0);
        
        
        if dislay_
            fprintf('  iteration %3d -> erreur = %.3e  \r',i,errorpf(i));
        end
        
        l0=l;U0=U;alpha0=alpha;
        l0tilde = ltilde;
        if errorpf(i)<paramradial.pfixtol
            break
        end
        
    end %%%%%%%%%%%%%%%% END BOUCLE POINT FIXE
    
    
    u = u + U.*l;
    if paramradial.adjoint
        utilde = utilde+Utilde*ltilde;
    end
    if ~paramradial.update && ~paramradial.updatefactor %%%%%%%%%%%%% NO UPDATE OF STOCHASTIC FUNCTIONS
        
        bu = bu - (A*U)*l;
        result.rayg{j}=result.rayg{j}(i);
        ee=result.raygini;
        for jj=1:j
            ee=ee+result.rayg{jj};
        end
        result.errorraygi(j)=sqrt(abs(result.rayg{j})/ee);
        
        radalpha=[radalpha, sparse(size(radalpha,1),1);...
            sparse(1,size(radalpha,2)),norm(l)];
        
    elseif paramradial.updatefactor
        V = double(getV(u));
        l = getL(u);
        
        fU = zeros(getm(u),1);
        aU = zeros(getm(u),getm(u));
        for kk=1:getm(u)
            fU(kk) = V(:,kk)'*expect(b,l(kk));
            for ll=1:1:getm(u)
                aU(kk,ll) = expect(V(:,kk)'*A*V(:,ll),l(kk),l(ll));
            end
        end
        fac = aU\fU;
        l= diag(fac)*l;
        result.Rayg{j} = double(full(expectmtimes(l,(V'*b)')));
        result.rayg{j} = trace(result.Rayg{j});
        u = PCRADIALMATRIX(V,[n,1],l);
        bu = b - A*u;
        if j==1
            result.errorrayg(j) = sqrt(abs(result.rayg{j}-result.raygini)/((result.rayg{j})));
        else
            result.errorrayg(j)=sqrt(abs(result.rayg{j}-result.rayg{j-1})/(result.rayg{j}));
        end
        
    else   %%%%%%%%%%%%%% UPDATING STOCHASTIC FUNCTIONS
        V = double(getV(u));
        fU = V'*b;
        aU = V'*A*V;
        
        [l,flag]=localstosolverupdate(aU,fU);
        l = expand(l);
        
        result.Rayg{j} = double(full(expectmtimes(l,fU')));
        result.rayg{j} = trace(result.Rayg{j});
        bsub = b;
        for kk=1:getparam(GSD,'subspaceiteration')
            if dislay_
                fprintf('subspace iteration %d \n',kk)
            end
            subspaceiteration
        end
        
        
        
        u = PCRADIALMATRIX(V,[n,1],l);
        bu = b - A*u;
        
        if j==1
            result.errorrayg(j) = sqrt(abs(result.rayg{j}-result.raygini)/((result.rayg{j})));
        else
            result.errorrayg(j)=sqrt(abs(result.rayg{j}-result.rayg{j-1})/(result.rayg{j}));
        end
        
    end %%%%%%%%%%%%% FIN UPDATING
    
    
    if getparam(GSD,'righthandSD');
        bu = righthandSD(GSD,bu);
    end
    
    
    % CALCUL ERREUR
    switch errorindicator
        case 'residual'
            result.error(j) = mynorm(bu)/normb;
        case 'rayleigh'
            if ~paramradial.update &&  ~paramradial.updatefactor
                result.error(j)= result.errorraygi(j)  ;
            else
                result.error(j)= result.errorrayg(j)  ;
            end
        case 'reference'
            result.error(j)=mynorm(u-uref)/mynorm(uref);
        case 'none'
            result.error(j)=1;
    end
    if ischarin('testerror',varargin)
        % result.errorA(j)=mynorm(u-uref,A)/mynorm(uref,A);
        result.errorres(j)=mynorm(bu)/normb;
    end
    
    
    if dislay_
        fprintf('  %d fonctions -> erreur = %.3e \n',result.nfonctions,result.error(j))
        if ischarin('testerror',varargin)
            fprintf('  %d fonctions -> erreurres = %.3e \n',result.nfonctions,result.errorres(j))
        end
    end
    
    result.addedfunctions = j;
    result.time(j)=etime(clock,clock0);
    
    if result.error(j)<paramradial.tol
        if strcmp(errorindicator,'rayleigh')
            u = normsort(u);
            u=truncate(u,length(u)-1);
            result.nfonctions = length(u);
            result.addedfunctions = j-1;
        end
        hasconv=true;
        fprintf('  Convergence avec %3d fonctions -> erreur = %.3e \n',result.nfonctions,result.error(j))
        break
    end
    
    hasconv=false;
    
end %%%%%%%%%%%%%%%% END BOUCLE NOUVELLES FONCTIONS


result.totaltime=etime(clock,clock0);



if getparam(GSD,'subspaceiteration')>0
    bsub = b;
    l = getL(u);
    for kk=1:getparam(GSD,'subspaceiteration')
        if dislay_
            fprintf('subspace iteration %d \n',kk)
        end
        
        subspaceiteration
        
        result.raygsub{j}(kk)=rayg;
        errorpf(kk)=abs((result.raygsub{j}(kk)-(kk>1)*result.raygsub{j}(kk-(kk>1)))/result.raygsub{j}(kk));
        
    end
    
    
    u = PCRADIALMATRIX(V,[n,1],l);
    bu = b - A*u;
    
    
    % CALCUL ERREUR
    switch errorindicator
        case 'residual'
            result.error(j) = mynorm(bu)/normb;
        case 'rayleigh'
            if ~paramradial.update
                result.error(j)= result.errorraygi(j)  ;
            else
                result.error(j)= result.errorrayg(j)  ;
            end
        case 'reference'
            
            result.error(j)=mynorm(u-uref)/mynorm(uref);
    end
    if ischarin('testerror',varargin)
        % result.errorA(j)=mynorm(u-uref,A)/mynorm(uref,A);
        result.errorres(j)=mynorm(bu)/normb;
    end
    
    if dislay_
        fprintf('  %d fonctions -> erreur = %.3e \n',result.nfonctions,result.error(j))
        if ischarin('testerror',varargin)
            fprintf('  %d fonctions -> erreurres = %.3e \n',result.nfonctions,result.errorres(j))
        end
    end
    
    
    
end



%% Finalisation : decomposition spectrale de la solution

if getparam(GSD,'finalunique')
    tolfactSD = getparam(GSD,'finaluniquefacttol');
    u=uniqueV(u,paramradial.tol*tolfactSD);
    
elseif getparam(GSD,'finalSD')
    tolfactSD = getparam(GSD,'finalSDfacttol');
    u = spectral_decomposition(u,'tol',paramradial.tol*tolfactSD,...
        'nbfoncmax',getm(u));
    
end



%% Recalcul de la base stochastique finale
if getparam(GSD,'finalupdate') %%%% FINAL UPDATING
    U = double(getV(u));
    fU = U'*b;
    aU = U'*A*U;
    [l,flag]=localstosolverupdate(aU,fU);
    
    u = PCRADIALMATRIX(U,[n,1],l);
    bu = b - A*u;
    switch errorindicator
        case 'reference'
            result.error(j)=mynorm(u-uref)/mynorm(uref);
        otherwise
            result.error(j) = mynorm(bu)/normb;
    end
    
    disp('  Reactualisation')
    fprintf('  nb fonctions %3d - erreur  : %3d \n',result.nfonctions,result.error(j))
end %%%% END FINAL UPDATING


%% Verdict sur la convergence
if ~hasconv && nbfoncmax>0
    fprintf('  Non-convergence avec %3d fonctions -> erreur = %.3e \n',result.nfonctions,result.error(j))
end

if dislay_
    fprintf('Elapsed time is = %.3f seconds', etime(clock,clock0))
end

