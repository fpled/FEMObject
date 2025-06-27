function [u,result] = solve_alterne_V2(A,b,solver,varargin)


% Initialisation param
if ( length(varargin)==1 && isa(varargin{1},'HSEPSOLVER') )
    % Dans ce cas, l'appel est realise depuis un HSM :
    param=solver;
else
    param = getparam(solver);
end

% Probleme en residu : stoquer (A/b) pour calcul d'erreur
if param.residual
    Apb=A;
    bpb=b;
    b=A'*b;
    A=A'*A;
else
    Apb=A;
    bpb=b;
end

dim = getdim(b);
n   = size(b);
bu  = b;
param.updatedim(n==1)=[]; % ???

% Initial guess
if ischarin('initialguess',varargin)
    u = getcharin('initialguess',varargin);
    bu = b - A*u;
    if param.righthandSD
        mm = getm(bu);
        if getdim(A)>2
            bu = multisvd(bu,'tol',param.tol/5,'maxorder',max(size(A)),'display',false);
            if param.display;fprintf('SVD of right-hand side : %d -> %d\n',mm,getm(bu));end
        else
            bu=svd(bu,param.tol/5);
            if param.display;fprintf('SVD of right-hand side : %d -> %d\n',mm,getm(bu)); end
        end
    end
else
    u = SEPMATRIX(dim);
end

% Initialisation
u0       = u;                        % utile pour || u_m+1 - u_m ||
utilde   = SEPMATRIX(dim);           % initialiser utilde
erriter  = zeros(1,param.maxiter);   % initialiser stockage
errorder = zeros(1,param.maxorder);  % initialiser stockage

% Metrique
if ~isempty(param.metric)
    Ametric = param.metric;
else
    Ametric = SEPMATRIX(cellfun(@(n)speye(n),num2cell(n),'UniformOutput',0));
end

% Initialiser normref
switch param.errorindicator
    case 'reference'
        try
            normref = normND(double(expand(param.reference)));
        catch
            normref = normND((param.reference));
        end
    case 'residual'
        normref = norm(bpb);
end

% Stoquer A sous un format convenable
AFR=mygathervector(A);
% Indicateur d'update
updated_a=1;
updated_d=1;

% Boucle Greedy
for i=1:param.maxorder
    
    % Affichage
    if param.display
        fprintf('%d',param.node);for cpt=1:param.depth, fprintf('      ');end
        fprintf('order #%d \n',i)
    end
    
    % Stoquer bu sous un format agreable
    if updated_d==1
        BFR=mygathervector(bu);
    else
        BFRnew=mygathervector(truncate(bu,bu.m-A.m+1 : bu.m));
        BFR=cellfun(@(b,bn) [b bn] , BFR , BFRnew ,'UniformOutput',0);
    end
    
    % Initialisation correction (U,alpha) / (Utilde,alphatilde)
    U0 = seprand(n);
    U = U0;
    Utilde = U0;
    alpha0 = 0;
    
    % Calcul de UAU et UB
    UAU = zeros(size(AFR{1},2),dim+1);
    UB  = zeros(size(BFR{1},2),dim+1);
    for d=1:dim
        tmp=Utilde.F{d} * U.F{1,d}';
        UAU(:,d)=tmp(:)'*AFR{d};
        UB(:,d)  = Utilde.F{d}' * BFR{d};
    end
    UAU(:,end) = A.alpha;
    UB(:,end)  = bu.alpha;
    
    % Alternating Least Square
    for kkk=1:param.maxiter
        
        % Boucle sur les dimensions
        for d=1:dim
            
            % Reduire bu
            db=1:dim+1;
            db(d)=[];
            ubd=prod(UB(:,db),2);
            S=BFR{d}*ubd;
            
            % Reduire A
            uaud=sparse(prod(UAU(:,db),2));
            M=AFR{d}*uaud;
            M=reshape(M,n(d),n(d));
            
            % Resoudre
            U.F{d} = M\S;
            
            % Normaliser
            alpha  = norm(U.F{d});
            U.F{d} = U.F{d}/alpha;
            
            % Actualiser Utilde
            if param.adjoint == 1
                % Peut etre optimise
                if ismember(d,param.adjointdim)
                    switch param.adjointtype
                        case 0
                            btilde = Ametric*(alpha*U);
                            Ubtilde = U'*btilde;
                            Ubtilde.F(:,d) = btilde.F(:,d);
                        case 1
                            btilde = Ametric*(alpha*U + u) ;
                            Ubtilde = U'*btilde;
                            Ubtilde.F(:,d) = btilde.F(:,d);
                        case 2
                            btilde = Ametric*(alpha*U + u) - (A')*utilde ;
                            Ubtilde = U'*btilde;
                            Ubtilde.F(:,d) = btilde.F(:,d);
                        case 3
                            btilde = Ametric*(alpha*U) - (A')*utilde ;
                            Ubtilde = U'*btilde;
                            Ubtilde.F(:,d) = btilde.F(:,d);
                            
                        otherwise
                            error('pas prevu')
                    end
                    
                    S = reduce(Ubtilde,d);
                    Utilde.F{d} = (M')\S;
                else
                    Utilde.F{d} = U.F{d};
                end
                alphatilde  = norm(Utilde.F{d});
                Utilde.F{d} = Utilde.F{d}/alphatilde;
            else
                alphatilde = alpha;
                Utilde     = U;
            end
            
            % Actualiser UAU/UB
            tmp=Utilde.F{d} * U.F{1,d}';
            UAU(:,d) = tmp(:)'*AFR{d};
            UB(:,d)  = Utilde.F{d}'*BFR{d};
            
        end % Fin boucle sur dimension
        
        % Critere en stagnation
        % A REVOIR : ajouter un carre ?
        erriter(kkk)=abs(alpha-alpha0)/(alpha+alpha0);
        alpha0=alpha;
        
        result.alpha{i}(kkk)=alpha;
        result.alphatilde{i}(kkk)=alphatilde;
        
        % Affichage
        if param.display && ~(param.node)
            fprintf('%d ',param.node);
            fprintf('iteration #%d - stagnation = %d\n',kkk,erriter(kkk))
        elseif param.display && param.node && kkk==1
            fprintf('%d',param.node);for cpt=1:param.depth, fprintf('      ');end
            fprintf('iteration #%d- %.1d ',kkk,erriter(kkk))
        elseif param.display && param.node && kkk~=1
            fprintf('#%d- %.1d ',kkk,erriter(kkk))
        end
        
        % Sortie si critere verifie
        if erriter(kkk)<param.itercrit
            break
        end
    end  % Fin ALS
    
    % Ajout correction
    u = u + alpha*U;
    if param.adjoint
        utilde = utilde + alphatilde*Utilde;
    end
    
    % Affichage
    if param.display && param.node && kkk~=1
        fprintf('\n')
    end
    
    %%%%%% Quelques update
    updated_a=0;
    updated_d=0;
    % update-dim
    if param.update>0 && getm(u)>1  && mod(i,param.updatestep)==0
        if param.adjoint
            [u,utilde]=solve_update(A,b,u,utilde,Ametric,param);
        else
            [u]=solve_update(A,b,u,[],Ametric,param);
        end
        bu = b - A*u;
        updated_a=1;
        updated_d=1;
    else
        bu = bu - alpha*A*U;
    end
    % update-tucker
    if param.updatetucker && getm(u)>1 && mod(i,param.updatestep)==0 %&& param.update==0
        if param.display; fprintf('%d ',param.node);for cpt=1:param.depth, fprintf('      ');end; fprintf(' updatetucker\n');end
        
        if param.adjoint && param.updateadjoint
            [u,utilde] = solve_updatetucker(A,b,u,utilde,param);
        else
            u = solve_updatetucker(A,b,u,[],param);
        end
        bu = b - A*u;
        updated_a=1;
        updated_d=1;
    end
    % update-alpha
    if param.alphaupdate && mod(i,param.updatestep)==0
        if param.display;fprintf('%d ',param.node);for cpt=1:param.depth, fprintf('      ');end;fprintf(' alphaupdate\n');end
        if param.adjoint && param.updateadjoint
            u = solve_alphaupdate(A,b,u,utilde);
        else
            u = solve_alphaupdate(A,b,u);
        end
        bu = b - A*u;
        updated_a=1;
        updated_d=0;
    end
    %%%%%% Quelques update FIN
    
    % Recompression du second membre
    if param.righthandSD && mod(i,param.righthandSDstep)==0
        mm = getm(bu);
        if getdim(A)>2
            bu = multisvd(bu,'tol',param.tol/5,'maxorder',max(size(A)),'display',false);
            if param.display;fprintf('%d ',param.node);for cpt=1:param.depth, fprintf('      ');end;fprintf('SVD of right-hand side : %d -> %d\n',mm,getm(bu));end
        else
            bu=svd(bu,param.tol/5);
            if param.display;fprintf('%d ',param.node);for cpt=1:param.depth, fprintf('      ');end;fprintf('SVD of right-hand side : %d -> %d\n',mm,getm(bu)); end
        end
        updated_a=1;
        updated_d=1;
    end
    
    % Stoquer u_m
    if param.storeiter==1
        result.u{i}=u;
    end
    
    % Indicateur d'erreur
    switch param.errorindicator
        case 'reference'
            try
                errorder(i)=normND(expand(u)-double(expand(param.reference)))/normref;
            catch
                errorder(i)=normND((u)-(param.reference))/normref;
            end
        case 'residual'
            errorder(i)=norm_residu(Apb,u,bpb,0)/normref;
        case 'stagnation'
            errorder(i)=norm(u-u0)/norm(2*u);
        otherwise
            errorder(i)= sqrt(min(abs(u.alpha))^2/sum(u.alpha.^2));
    end
    if param.display
        fprintf('%d',param.node);for cpt=1:param.depth, fprintf('      ');end
        fprintf('order #%d - error = %d \n',i,errorder(i))
    end
    u0=u;
    
    % Sortie si convergence
    if errorder(i)<param.tol
        hasconv = 1;
        break
    end
end

% Affichage final
if (nargout==1)&&param.display;
    if errorder(i)>=param.tol
        fprintf('%d',param.node);for cpt=1:param.depth, fprintf('      ');end
        fprintf('---> no convergence - ')
    else
        fprintf('%d ',param.node);for cpt=1:param.depth, fprintf('      ');end
        fprintf('---> convergence    - ')
    end
    fprintf('order #%d - error = %d\n',i,errorder(i))
end

result.error = errorder;
% Stoquer utilde
if param.adjoint
    result.utilde = utilde;
end


end



function u=mygathervector(u)
dim=u.dim;
u=cellfun(@(a) a(:), u.F,'UniformOutput',0);
for d=1:dim
    u{1,d}=horzcat(u{:,d});
end
u=u(1,:);
end



















