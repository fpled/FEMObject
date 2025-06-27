function [u,result] = dsolve_power(GSD,T,b,A,B,u0,varargin)
% function u = dsolve_power(GSD,A,b,v,varargin)
% fonction solve GSD : resolution de Au=b
% A et b : double ou TIMEMATRIX
% v : pour la reutilisation (MULTIMATRIX ou TIMERADIALMATRIX ou double)

if nargin>=6 && ~isempty(u0)
    error('non prise en compte de la condition initiale')
end

display_ = getparam(GSD,'display');

fprintf('GSD solveur ... ')
if display_
    fprintf('\n')
end

paramradial.nbfoncmax = getparam(GSD,'nbfoncmax');
paramradial.tol = getparam(GSD,'tol');
paramradial.pfixtol = getparam(GSD,'pfixtol');
paramradial.pfixmax = getparam(GSD,'pfixmax');
paramradial.reuse = getparam(GSD,'reuse');
paramradial.update = getparam(GSD,'update');


errorindicator = getparam(GSD,'errorindicator');
uref = getcharin('reference',varargin);

if ~isempty(uref)
    errorindicator='reference';
end

[b,A,B] = init_resolution(T,b,A,B);
[b,A,B,PC]=pcsystemupdate(b,A,B,varargin{:});

n = size(A,1);
nt = length(T);
P=length(PC);

if getparam(GSD,'righthandSD');
    tolfactSD = getparam(GSD,'finalSDfacttol');
    
    if isa(b,'PCRADIALMATRIX')
        mm=length(b);
    else
        mm=100;
    end
    b = spectral_decomposition(b,'tol',paramradial.tol*tolfactSD,...
        'nbfoncmax',mm);
end


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



clock0 = clock;

%% D�finition des normes
Am = expect(A);
Am = 1/2*(Am+Am');
Am = speye(n);


if ~strcmp(errorindicator,'rayleigh')
    mynorm=@(b) norm(b);
    normb = mynorm(b);
end

%% D�finition des solveurs
toliter = getparam(GSD,'toliter');
directsolve = getparam(GSD,'direct');

if directsolve
    localstosolver = @(aU,fU) solve(aU,fU);
else
    localstosolver = @(aU,fU) cgs(aU,fU,toliter,[],'noupdate');
end


u = PCTIMEMATRIX(PCRADIALMATRIX([n,nt],PC),T,[n,1]);
bu = b;
result.errorini=1;

result.nfonctions = getm(u);
result.addedfunctions = 0;


radalpha = [];
V = cell(1,0);
DV = cell(1,0);


%% Construction des couples de fonctions
nbfoncmax = paramradial.nbfoncmax;

for j=1:nbfoncmax %%%%%%%%%%%%%%%% BOUCLE NOUVELLES FONCTIONS
    
    result.nfonctions = result.nfonctions + 1  ;
    
    % l0 = one(PC);
    l0 = ones(length(PC),1);l0=PCMATRIX(l0,[1,1],PC);
    alpha0=norm(l0);l0=l0/alpha0;
    
    for i=1:paramradial.pfixmax %%%%%%%%%%%%%%%% BOUCLE POINT FIXE
        
        f0=expectmtimes(bu,l0);
        A0 = expectmtimes(A,l0,l0);
        B0 = expectmtimes(B,l0,l0);
        U = dsolve(T,f0,A0,B0);
        
        % rayg(i)=U'*f0;
        alpha = norm(U);
        U = U /alpha ;
        DU = diff(U);
        
        fU = integratemtimes(U',bu);
        aU = integratemtimes(U',A,DU);
        bU = integratemtimes(U',B,U);
        aU = aU+bU;
        
        [l,flag]=localstosolver(aU,fU);
        
        
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
    
    V = [V,{U}];
    DV = [DV,{DU}];
    
    if ~paramradial.update %%%%%%%%%%%%% NO UPDATE OF STOCHASTIC FUNCTIONS
        u = u + U*(alpha*l);
        bu = bu - (A*DU+B*U)*(alpha*l);
        result.rayg{j}=result.rayg{j}(i);
        ee=0;
        for jj=1:j
            ee=ee+result.rayg{jj};
        end
        result.errorraygi(j)=sqrt(abs(result.rayg{j})/ee);
        
        radalpha=[radalpha, sparse(size(radalpha,1),1);...
            sparse(1,size(radalpha,2)),alpha];
        
        
    else   %%%%%%%%%%%%%% UPDATING STOCHASTIC FUNCTIONS
        clear aU bU fU;
        
        for ii=1:length(V)
            fU{ii} = integratemtimes(V{ii}',b);
            for jj=1:length(V)
                aU{ii}{jj} = integratemtimes(V{ii}',A,DV{jj});
                bU{ii}{jj} = integratemtimes(V{ii}',B,V{jj});
                aU{ii}{jj} = aU{ii}{jj}+bU{ii}{jj};
            end
            aU{ii} = horzcat(aU{ii}{:});
        end
        fU = vertcat(fU{:});
        aU = vertcat(aU{:});
        
        [l,flag]=localstosolver(aU,fU);
        
        bu= b;
        u = PCTIMEMATRIX(PCRADIALMATRIX([n,nt],PC),T,[n,1]);
        for ii=1:length(V)
            Vi = PCRADIALMATRIX(getvalue(V{ii}),[n,nt],l(ii));
            u= u + PCTIMEMATRIX(Vi,T,[n,1]);
            bu = bu - (A*DV{ii}+B*V{ii})*l(ii);
        end
        
        
        result.Rayg{j} = double(full(expect(fU,l)));
        result.rayg{j} = trace(result.Rayg{j});
        
        
        if j==1
            result.errorrayg(j) = sqrt(abs(result.rayg{j})/((result.rayg{j})));
        else
            result.errorrayg(j)=sqrt(abs(result.rayg{j}-result.rayg{j-1})/(result.rayg{j}));
        end
        
    end %%%%%%%%%%%%% FIN UPDATING
    
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
    
    
    if display_
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
            u = truncate(u,length(u)-1);
            result.nfonctions = length(u);
            result.addedfunctions = j-1;
        end
        hasconv=true;
        fprintf('  Convergence avec %3d fonctions -> erreur = %.3e \n',result.nfonctions,result.error(j))
        break
    end
    
    hasconv=false;
    
end %%%%%%%%%%%%%%%% END BOUCLE NOUVELLES FONCTIONS


%% Finalisation : decomposition spectrale de la solution

if getparam(GSD,'finalunique')
    tolfactSD = getparam(GSD,'finaluniquefacttol');
    u=uniqueV(u,paramradial.tol*tolfactSD);
    
elseif getparam(GSD,'finalSD')
    tolfactSD = getparam(GSD,'finalSDfacttol');
    u = spectral_decomposition(u,'tol',paramradial.tol*tolfactSD,...
        'nbfoncmax',length(u));
end



%% Verdict sur la convergence
if ~hasconv
    fprintf('  Non-convergence avec %3d fonctions -> erreur = %.3e \n',result.nfonctions,result.error(j))
end

if display_
    fprintf('Elapsed time is = %.3f seconds', etime(clock,clock0))
end
