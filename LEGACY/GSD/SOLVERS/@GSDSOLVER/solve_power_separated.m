function [u,result] = solve_power_separated(GSD,A,b,ureuse,varargin)
% function u = solve(GSD,A,b,v,varargin)
% fonction solve GSD : resolution de Au=b
% A et b : PCMATRIX ou PCRADIALMATRIX
% v : pour la reutilisation (MULTIMATRIX ou PCRADIALMATRIX ou double)

display_ = getparam(GSD,'display');

fprintf('\n GSD solveur ... ')
if display_
    fprintf('\n')
end

paramradial.nbfoncmax = getparam(GSD,'nbfoncmax');
paramradial.tol = getparam(GSD,'tol');
paramradial.pfixtol = getparam(GSD,'pfixtol');
paramradial.pfixmax = getparam(GSD,'pfixmax');
paramradial.reuse = getparam(GSD,'reuse');
paramradial.update = getparam(GSD,'update');
paramradial.orthoduringpfix = getparam(GSD,'orthoduringpfix');

errorindicator = getparam(GSD,'errorindicator');
uref = getcharin('reference',varargin);

if ~isempty(uref)
    errorindicator='reference';
end

n=size(A,1);

errorpf=zeros(1,paramradial.pfixmax);
result.nfonctions = 0 ;
result.addedfunctions = 0;

PC = getPC(A);
l0 = zeros(1,PC);
l = zeros(1,PC);
U = sparse(n,1);
stoU = sparse(n,0);
stol = cell(1,0);
U0 = sparse(n,1);
radalpha = [];


%% D�finition des normes
Am = speye(n);
if ~strcmp(errorindicator,'rayleigh')
    mynorm=@(b) norm(b);
    normb = mynorm(b);
end

%% D�finition des solveurs
toliter = getparam(GSD,'toliter');
tolupdate = getparam(GSD,'tolupdate');

localstosolver = getparam(GSD,'localstosolver');

localstosolveriter = @(A,b) solve(localstosolver,A,b,toliter);
localstosolverupdate = @(A,b)  solve(localstosolver,A,b,tolupdate);

clock0 = clock;

result.totaltime=etime(clock,clock0);


u = PCTPRADIALMATRIX([n,1],PC);
result.raygini=0;
result.errorini=1;


%% Construction des couples de fonctions
nbfoncmax = paramradial.nbfoncmax;

for j=1:nbfoncmax %%%%%%%%%%%%%%%% BOUCLE NOUVELLES FONCTIONS
    
    result.nfonctions = result.nfonctions + 1  ;
    
    switch getparam(GSD,'inittype')
        case 'one'
            l0 = one(PC);
        case 'random'
            l0 = rand(PC);
        otherwise
            error('pas defini')
            
    end
    
    alpha0=norm(l0);l0=l0/alpha0;
    
    for i=1:paramradial.pfixmax %%%%%%%%%%%%%%%% BOUCLE POINT FIXE
        
        f0=expect(b,l0);
        for kk=1:getm(u)
            f0 = f0 - expect(A*getV(u,kk),l0,getL(u,kk));
        end
        L0 = expect(A,l0,l0);
        U = L0\f0;
        
        if paramradial.orthoduringpfix
            lllll
            V = double(getV(u));
            if size(V,2)>=1
                U = U - V*(V'*Am*U);
            end
        end
        alpha = sqrt(U'*Am*U);
        U = U /alpha ;
        fU = U'*b;
        for kk=1:getm(u)
            fU = fU - (U'*A*getV(u,kk))*getL(u,kk);
        end
        
        
        aU = U'*A*U;
        
        
        [l,flag]=localstosolveriter(aU,fU);
        
        result.rayg{j}(i) = full(expectmtimes(l,fU'));
        
        errorpf(i)=abs((result.rayg{j}(i)-(i>1)*result.rayg{j}(i-(i>1)))/result.rayg{j}(i));
        
        if display_
            fprintf(' GSD iteration %3d -> erreur = %.3e  \r',i,errorpf(i));
        end
        
        l0=l;U0=U;alpha0=alpha;
        
        if errorpf(i)<paramradial.pfixtol
            break
        end
        
    end %%%%%%%%%%%%%%%% END BOUCLE POINT FIXE
    
    u = u + PCTPRADIALMATRIX(U,size(U),{l});
    
    if ~paramradial.update  %%%%%%%%%%%%% NO UPDATE OF STOCHASTIC FUNCTIONS
        
        result.rayg{j}=result.rayg{j}(i);
        ee=result.raygini;
        for jj=1:j
            ee=ee+result.rayg{jj};
        end
        result.errorraygi(j)=sqrt(abs(result.rayg{j})/ee);
        
        radalpha=[radalpha, sparse(size(radalpha,1),1);...
            sparse(1,size(radalpha,2)),norm(l)];
        
        
    else   %%%%%%%%%%%%%% UPDATING STOCHASTIC FUNCTIONS
        
        V = double(getV(u));
        fU = V'*b;
        aU = V'*A*V;
        
        if display_
            fprintf('Actualisation des variables aleatoires\n')
        end
        [l,flag]=localstosolverupdate(aU,fU);
        
        result.Rayg{j} = double(full(expectmtimes(l,fU')));
        result.rayg{j} = trace(result.Rayg{j});
        
        u = PCTPRADIALMATRIX(V,[n,1],l);
        
        if j==1
            result.errorrayg(j) = sqrt(abs(result.rayg{j}-result.raygini)/((result.rayg{j})));
        else
            result.errorrayg(j)=sqrt(abs(result.rayg{j}-result.rayg{j-1})/(result.rayg{j}));
        end
        
    end %%%%%%%%%%%%% FIN UPDATING
    
    
    
    % CALCUL ERREUR
    switch errorindicator
        case 'residual'
            result.error(j) = mynorm(b-A*u)/normb;
        case 'rayleigh'
            if ~paramradial.update
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
    
    
    if display_
        fprintf(' GSD %d fonctions -> erreur = %.3e \n',result.nfonctions,result.error(j))
        if ischarin('testerror',varargin)
            fprintf(' GSD %d fonctions -> erreurres = %.3e \n',result.nfonctions,result.errorres(j))
        end
    end
    
    result.addedfunctions = j;
    result.time(j)=etime(clock,clock0);
    
    if result.error(j)<paramradial.tol
        hasconv=true;
        fprintf(' GSD Convergence avec %3d fonctions -> erreur = %.3e \n',result.nfonctions,result.error(j))
        break
    end
    
    hasconv=false;
    
end %%%%%%%%%%%%%%%% END BOUCLE NOUVELLES FONCTIONS


result.totaltime=etime(clock,clock0);


%% Verdict sur la convergence
if ~hasconv && nbfoncmax>0
    fprintf('  Non-convergence avec %3d fonctions -> erreur = %.3e \n',result.nfonctions,result.error(j))
end

if display_
    fprintf('Elapsed time is = %.3f seconds', etime(clock,clock0))
end
