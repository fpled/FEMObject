function [u,result] = funsolve(GSD,Aexpect,AU,b,ureuse,varargin)
% function u = solve(GSD,A,b,v,varargin)
% fonction solve GSD : resolution de Au=b
% A et b : PCMATRIX ou PCRADIALMATRIX
% v : pour la reutilisation (MULTIMATRIX ou PCRADIALMATRIX ou double)

fprintf('\nGSD solveur ... \n')
uref = getcharin('reference',varargin);
paramradial.nbfoncmax = getparam(GSD,'nbfoncmax');
paramradial.tol = getparam(GSD,'tol');
paramradial.pfixtol = getparam(GSD,'pfixtol');
paramradial.pfixmax = getparam(GSD,'pfixmax');
paramradial.reuse = getparam(GSD,'reuse');
errorindicator = getparam(GSD,'errorindicator');
update = getparam(GSD,'update');

if ~isempty(uref)
    errorindicator='reference';
end

tolpcg = getparam(GSD,'reuse');
if paramradial.reuse
    qreuse = getcharin('reuse',varargin);
    varargin = delcharin('reuse',varargin);
end
tolpcg = getparam(GSD,'toliter');
directsolve = getparam(GSD,'direct');
display_ = getparam(GSD,'display');

[Aexpect,b,PC]=pcsystemupdate(Aexpect,b,varargin{:});
if ~isa(AU,'inline') & ~isa(AU,'function_handle')
    calc_AU = @(U) AU*U;
else
    calc_AU = AU;
end
P = getP(PC);
n=size(b,1);

errorpf=zeros(1,paramradial.pfixmax);
result.error=zeros(1,paramradial.nbfoncmax);
result.errorres=zeros(1,paramradial.nbfoncmax);
result.errorrayg=zeros(1,paramradial.nbfoncmax);
result.errorraygi=zeros(1,paramradial.nbfoncmax);
result.time = zeros(1,paramradial.nbfoncmax);
result.rayg=cell(1,paramradial.nbfoncmax);
result.Rayg=cell(1,paramradial.nbfoncmax);
result.nfonctions = 0 ;
result.addedfunctions = 0;

l0 = zeros(1,PC);
l = zeros(1,PC);
U = sparse(n,1);
stoU = sparse(n,0);
U0 = sparse(n,1);
radalpha = [];

Am = Aexpect([],[]);
Am = 1/2*(Am+Am');

tic

time0 = toc;

if paramradial.reuse
    if isa(ureuse,'PCRADIALMATRIX')
        if length(ureuse)==0
            paramradial.reuse=0;
        else
            ureuse = getV(ureuse);
        end
    elseif ~isa(ureuse,'MULTIMATRIX')
        if isa(ureuse,'PCMATRIX')
            fprintf('  on ne peut reutiliser une PCMATRIX -> utiliser PCRADIALMATRIX ou MULTIMATRIX\n')
            paramradial.reuse=0;
        end
        if isa(ureuse,'double') & normest(ureuse)<eps
            fprintf('  vecteur de norme negligeable -> pas d''initialisation\n')
            paramradial.reuse=0;
        end
    end
end

if ~paramradial.reuse    % PAS D'INITIALISATION
    u = PCRADIALMATRIX([n,1],PC);
    bini = b;
    
else  % PHASE D'INITIALISATION
    
    U = double(ureuse);
    fU = U'*b;
    aU = U'*calc_AU(U);
    if directsolve
        l = solve(aU,fU);
    else
        [l,flag] = cgs(aU,fU,tolpcg,[],'noupdate');
    end
    
    u = PCRADIALMATRIX(U,[n,1],l);
    bini = b - calc_AU(u);
    
    switch errorindicator
        case 'reference'
            result.errorini=norm(u-uref)/norm(uref);
        otherwise
            result.errorini=norm(bini)/norm(b);
    end
    result.nfonctions = length(u);
    fprintf('  Initialisation : %3d fonctions - erreur  : %3d \n',result.nfonctions,result.errorini)
    if result.errorini<getparam(GSD,'tolini')
        disp('  Convergence apres initialisation')
        return
    end
    
end

bu = bini;

nbfoncmax = paramradial.nbfoncmax;

for j=1:nbfoncmax %%%%%%%%%%%%%%%% BOUCLE NOUVELLES FONCTIONS
    
    result.nfonctions = result.nfonctions + 1  ;
    
    % l0 = one(PC);
    l0 = ones(length(PC),1);l0=PCMATRIX(l0,[1,1],PC);
    alpha0=norm(l0);l0=l0/alpha0;
    
    for i=1:paramradial.pfixmax %%%%%%%%%%%%%%%% BOUCLE POINT FIXE
        
        f0=expect(bu,l0);
        L0 = Aexpect(l0,l0);
        U = L0\f0;
        % rayg(i)=U'*f0;
        alpha = sqrt(U'*Am*U);
        U = U /alpha ;
        fU = U'*bu;
        aU = U'*calc_AU(U);
        % if directsolve
        l = solve(aU,fU);
        % else
        % [l,flag] = cgs(aU,fU,tolpcg,[],'noupdate');
        % end
        
        result.rayg{j}(i) = expect(l,PCMATRIX(fU));
        alpha = norm(l);
        l=l/alpha;
        % errorpf(i)=double(sqrt(...
        %   (alpha0^2+alpha^2-2*alpha0*alpha*(U0'*Am*U)*(double(l0)*double(l)'))/...
        %                 (alpha^2)));
        errorpf(i)=abs((result.rayg{j}(i)-(i>1)*result.rayg{j}(i-(i>1)))/result.rayg{j}(i));
        
        if display_
            fprintf('  iteration %3d -> erreur = %.3e  \r',i,errorpf(i));
        end
        
        l0=l;U0=U;alpha0=alpha;
        
        if errorpf(i)<paramradial.pfixtol
            break
        end
        
    end %%%%%%%%%%%%%%%% END BOUCLE POINT FIXE
    
    u = u + U.*(alpha*l);
    
    if ~update  %%%%%%%%%%%%% NO UPDATE OF STOCHASTIC FUNCTIONS
        bu = bu - (calc_AU(U))*(alpha*l);
        ee=0;
        for jj=1:j
            ee=ee+result.rayg{jj}(end);
        end
        result.errorraygi(j)=abs(result.rayg{j}(end))/ee;
        
        radalpha=[radalpha, sparse(size(radalpha,1),1);...
            sparse(1,size(radalpha,2)),alpha];
        
        
    else   %%%%%%%%%%%%%% UPDATING STOCHASTIC FUNCTIONS
        V = double(getV(u));
        fU = V'*b;
        aU = V'*calc_AU(V);
        
        if directsolve
            l = solve(aU,fU);
        else
            [l,flag] = cgs(aU,fU,tolpcg,[],'noupdate');
        end
        
        result.Rayg{j} = double(full(expect(fU,l)));
        
        u = PCRADIALMATRIX(V,[n,1],l);
        bu = b - calc_AU(u);
        
        result.errorrayg(j)=abs(trace(result.Rayg{j})-(j>1)*trace(result.Rayg{j-(j>1)}))/trace(result.Rayg{j});
        
    end %%%%%%%%%%%%% FIN UPDATING
    
    % CALCUL ERREUR
    switch errorindicator
        case 'residual'
            result.error(j) = norm(bu)/norm(b);
        case 'rayleigh'
            if ~update
                result.error(j)= result.errorraygi(j)  ;
            else
                result.error(j)= result.errorrayg(j)  ;
            end
        case 'reference'
            result.error(j)=norm(u-uref)/norm(uref);
    end
    if ischarin('testerror',varargin)
        % result.errorA(j)=norm(u-uref,A)/norm(uref,A);
        result.errorres(j)=norm(bu)/norm(b);
    end
    
    
    if display_
        fprintf('  %d fonctions -> erreur = %.3e \n',result.nfonctions,result.error(j))
        if ischarin('testerror',varargin)
            fprintf('  %d fonctions -> erreurres = %.3e \n',result.nfonctions,result.errorres(j))
        end
    end
    
    result.addedfunctions = j;
    result.time(j)=toc-time0;
    if result.error(j)<paramradial.tol
        fprintf('  Convergence avec %3d fonctions -> erreur = %.3e \n',result.nfonctions,result.error(j))
        return
    end
    
end %%%%%%%%%%%%%%%% END BOUCLE NOUVELLES FONCTIONS


if getparam(GSD,'finalupdate') %%%% FINAL UPDATING
    U = double(getV(u));
    fU = U'*b;
    aU = U'*calc_AU(U);
    if directsolve
        l = solve(aU,fU);
    else
        [l,flag] = cgs(aU,fU,tolpcg,[],'noupdate');
    end
    
    u = PCRADIALMATRIX(U,[n,1],l);
    bu = b - calc_AU(u);
    switch errorindicator
        case 'reference'
            result.error(j)=norm(u-uref)/norm(uref);
        otherwise
            result.error(j) = norm(bu)/norm(b);
    end
    
    disp('  Reactualisation')
    fprintf('  nb fonctions %3d - erreur  : %3d \n',result.nfonctions,result.error(j))
end %%%% END FINAL UPDATING

fprintf('  Non-convergence avec %3d fonctions -> erreur = %.3e \n',result.nfonctions,result.error(j))


toc