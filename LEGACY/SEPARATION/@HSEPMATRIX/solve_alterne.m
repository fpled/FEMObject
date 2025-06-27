function [u,result] = solve_alterne(A,b,solver,varargin)
% function [u,result] = solve_alterne(A,b,solver,varargin)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~(length(varargin)==1 && isa(varargin{1},'PARAMETERS') ) %Premier passage
    param  = solver(1);
else % Recurence...
    param  = solver;
    solver = varargin{1};
end
Tp=param.tree;
SubNodes=find(getconnect(Tp)==param.node);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ta  = A.tree;
dim = A.dim;
n   = size(A);
bu  = b;
repone = find(n==1);         % ??
param.updatedim(repone)=[];  % ??

if ischarin('initialguess',varargin)
    error('Non implemente');
    %     u = getcharin('initialguess',varargin);
    %     bu = b - A*u;
    %     if param.righthandSD
    %         mm = getm(bu);
    %         if getdim(A)>2
    %             bu = multisvd(bu,'tol',param.tol/5,'maxorder',max(size(A)),'display',false);
    %             if param.display;fprintf('SVD of right-hand side : %d -> %d\n',mm,getm(bu));end
    %         else
    %             bu=svd(bu,param.tol/5);
    %             if param.display;fprintf('SVD of right-hand side : %d -> %d\n',mm,getm(bu)); end
    %         end
    %     end
else
    u = HSEPMATRIX(Ta);
end

u0=u;
utilde = HSEPMATRIX(Ta);
hasconv = 0; % ??
erriter  = zeros(1,param.maxiter);
errorder = zeros(1,param.maxorder);

if ~isempty(param.metric)
    Ametric = param.metric;
else
    Ametric = heye(Ta,n);
end


switch param.errorindicator
    case 'reference'
        error('Non implemente');
        %         try
        %             normref = normND(double(expand(param.reference)),Ametric);
        %         catch
        %             normref = normND((param.reference),Ametric);
        %         end
    case 'residual'
        normref = norm(b);
end

for i=1:param.maxorder
    
    if param.display
        fprintf('%d ',param.node);for cpt=1:param.depth, fprintf('      ');end
        fprintf('order #%d \n',i)
    end
    
    % Initialisation de U, la correction :
    switch param.inittype
        case 'rand'
            U0 = hseprand(Ta,n);
        case 'one'
            U0 = hsepone(Ta,n);
    end
    U0tilde = U0;
    alpha0=0;
    U = U0;
    Utilde = U0tilde;
    
    % Calcul rapide de Utilde'*A*U
    AU = A*U;
    UAUF = cellfun(@fastprodscal,Utilde.F(ones(AU.m,1),:),AU.F);
    UAUA = Utilde.alpha(1)*AU.alpha;
    
    % Calcul rapide de Utilde'*(b-U)
    UbuF = cellfun(@fastprodscal,Utilde.F(ones(bu.m,1),:),bu.F);
    UbuA = Utilde.alpha(1)*bu.alpha;
    
    for kkk=1:param.maxiter    %%%%  ITERATION  ALTERNEE
        for j=1:dim %%%% BOUCLE DIMENSION j
            M = reduce(UAUF,UAUA,j,A.F(:,j));
            S = reduce(UbuF,UbuA,j,bu.F(:,j));
            U.F{j}    = solve_alterne(M,S,...
                solver(SubNodes(j)),...
                setparam(solver,'node',SubNodes(j)));
            alpha     = norm(U.F{j});
            U.F{j}    = U.F{j}*(1/alpha);
            if param.adjoint == 1
                % Non implemente
            else
                alphatilde = alpha;
                Utilde = U;
            end
            % Actualisation operateur
            AU.F(:,j) = cellfun(@mtimes,A.F(:,j),U.F(ones(AU.m,1),j),'UniformOutput',0);
            UAUF(:,j) = cellfun(@fastprodscal,Utilde.F(ones(AU.m,1),j),AU.F(:,j));
            UAUA      = Utilde.alpha(1)*AU.alpha;
            
            % Actualisation 2nd membre
            UbuF(:,j) = cellfun(@fastprodscal,Utilde.F(ones(bu.m,1),j),bu.F(:,j));
            UbuA      = Utilde.alpha(1)*bu.alpha;
        end %%%% END BOUCLE DIMENSION j
        
        % Calcul de l'erreur en stagnation :
        erriter(kkk)=abs(alpha-alpha0)/(alpha+alpha0);
        result.alpha{i}(kkk)=alpha;
        result.alphatilde{i}(kkk)=alphatilde;
        if param.display
            fprintf('%d ',param.node);for cpt=1:param.depth, fprintf('      ');end
            fprintf('iteration #%d - stagnation = %d\n',kkk,erriter(kkk))
        end
        % Detection de stabilisation :
        if erriter(kkk)<param.itercrit
            break
        end
        
        % Chaise musicale : ancien:=nouveau
        alpha0=alpha;
        alpha0tilde=alphatilde;
        U0 = U;
        U0tilde=Utilde;
        
    end  %%%%  END ITERATION ALTERNEE
    
    if param.display
        fprintf('%d ',param.node);for cpt=1:param.depth, fprintf('      ');end
        for dd=1:U.dim
            fprintf('Dim %d Rang %d     ',dd,U.F{1,dd}.m) ;
        end
        fprintf('\n');
    end
    
    % AJOUT DU NOUVEAU RANG :
    u = u + alpha*U;
    
    % Si adjoint
    if param.adjoint
        utilde = utilde + alphatilde*Utilde;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%% UPDATE %%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % UpdateDim
    if param.update>0 && i>1 && mod(i,param.updatestep)==0
    else
        bu = bu - alpha*AU;
    end
    % UpdateDim
    if param.alphaupdate && mod(i,param.updatestep)==0
        if param.display;fprintf('%d ',param.node);for cpt=1:param.depth, fprintf('      ');end;fprintf(' alphaupdate\n');end
        if param.adjoint && param.updateadjoint
            u = solve_alphaupdate(A,b,u,utilde);
        else
            try
                u = solve_alphaupdate(A,b,u);
            catch
                caca=1;
            end
        end
        bu = b - A*u;
    end
    % Tucker...
    if param.updatetucker && dim==2 && mod(i,param.updatestep)==0 && param.update==0
    end
    
    if param.righthandSD && mod(i,param.righthandSDstep)==0
    end
    
    
    % Calcul de l'erreur de l'approximation (au rang i)
    switch param.errorindicator
        case 'reference'
            try
                errorder(i)=normND(expand(u)-double(expand(param.reference)))/normref;
            catch
                errorder(i)=normND((u)-(param.reference))/normref;
            end
        case 'residual'
            errorder(i)=norm(bu)/normref;
        case 'stagnation'
            errorder(i)=norm(u-u0)/normND(u+u0);
        otherwise
            errorder(i)= sqrt(min(abs(u.alpha))^2/sum(u.alpha.^2));
    end
    
    
    if param.node==1%root==1 || param.display
        fprintf('%d ',param.node);for cpt=1:param.depth, fprintf('      ');end
        fprintf('order #%d - error = %d \n',i,errorder(i))
    end
    u0=u;
    
    if param.adjoint
        result.utilde = utilde;
    end
    
    if errorder(i)<param.tol
        hasconv = 1;
        break
    end
    
end

if (nargout==1)&&param.display;
    if errorder(i)>=param.tol
        fprintf('%d ',param.node);for cpt=1:param.depth, fprintf('      ');end
        fprintf('---> no convergence - ')
    else
        fprintf('%d ',param.node);for cpt=1:param.depth, fprintf('      ');end
        fprintf('---> convergence    - ')
    end
    fprintf('order #%d - error = %d\n',i,errorder(i))
end
result.error = errorder;
if param.adjoint
    result.utilde = utilde;
end


function M = reduce(F,A,j,Z)
% Compression de l'info :
%   produit des autres dimension
%   somme sur les rangs
F(:,j)=1;
a = (A).*prod(F,2)';
M = a(1)*Z{1};
for i=2:length(a)
    M = M + a(i)*Z{i};
end
return


