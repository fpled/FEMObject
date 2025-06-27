function [u,result] = solve_alterne_non_lineaire(A,b,solver,varargin)
% function [u,result] = solve_alterne_non_lineaire(A,b,solver,varargin)

dim = getdim(A);
if ( length(varargin)==1 && isa(varargin{1},'HSEPSOLVER') )
    % Dans ce cas, l'appel est realise depuis un HSM :
    param=solver;
else
    param = getparam(solver);
end

if param.residual
    solver = setparam(solver,'residual',false);
    [u,result] = solve_alterne(A'*A,A'*b,solver,varargin{:});
    return
end

n = size(b);
repone = find(n==1);
param.updatedim(repone)=[];
bu = b;

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

u0=u;
utilde = SEPMATRIX(dim);
hasconv = 0;
erriter = zeros(1,param.maxiter);
errorder = zeros(1,param.maxorder);

if ~isempty(param.metric)
    Ametric = param.metric;
else
    Ametric = cell(1,dim);
    for j=1:dim
        Ametric{j} = speye(n(j));
    end
    Ametric = SEPMATRIX(Ametric);
end

% if dim==2
%     param.maxorder = min(min(n),param.maxorder);
% end

switch param.errorindicator
    case 'reference'
        try
            normref = normND(double(expand(param.reference)));
        catch
            normref = normND((param.reference));
        end
    case 'residual'
        normref = norm(b);
end

for i=1:param.maxorder
    
    
    if param.display
        fprintf('%d',param.node);for cpt=1:param.depth, fprintf('      ');end
        fprintf('order #%d \n',i)
    end
    
    switch param.inittype
        case 'rand'
            U0 = seprand(n);
        case 'one'
            U0 = sepone(n);
    end
    
    U0tilde = U0;
    alpha0=0;
    
    U = U0;
    Utilde = U0tilde;
    AU = A*U;
    UAU = Utilde'*A*U;
    Ubu = Utilde'*bu;
    
    for kkk=1:param.maxiter    %%%%  ITERATION  ALTERNEE
        
        for j=1:dim %%%% BOUCLE DIMENSION j
            
            dNewtonErr = param.bNlSolver + param.bNlSolver;
            Up = U;
            while (dNewtonErr>=param.bNlSolver && dNewtonErr>=0)
                UAU.F(:,j)= A.F(:,j);
                M = reduce(UAU,j);
                Ubu.F(:,j) = bu.F(:,j);
                S = reduce(Ubu,j);
                
                U.F{j} = M\S;
                alpha = norm(U.F{j});
                
                U.F{j} = U.F{j}/alpha;
                
                if param.adjoint == 1
                    if ismember(j,param.adjointdim)
                        switch param.adjointtype
                            case 0
                                btilde = Ametric*(alpha*U);
                                Ubtilde = U'*btilde;
                                Ubtilde.F(:,j) = btilde.F(:,j);
                            case 1
                                btilde = Ametric*(alpha*U + u) ;
                                Ubtilde = U'*btilde;
                                Ubtilde.F(:,j) = btilde.F(:,j);
                            case 2
                                btilde = Ametric*(alpha*U + u) - (A')*utilde ;
                                Ubtilde = U'*btilde;
                                Ubtilde.F(:,j) = btilde.F(:,j);
                            case 3
                                btilde = Ametric*(alpha*U) - (A')*utilde ;
                                Ubtilde = U'*btilde;
                                Ubtilde.F(:,j) = btilde.F(:,j);
                                
                            otherwise
                                error('pas prevu')
                        end
                        
                        S = reduce(Ubtilde,j);
                        Utilde.F{j} = (M')\S;
                    else
                        Utilde.F{j} = U.F{j};
                    end
                    alphatilde = norm(Utilde.F{j});
                    Utilde.F{j} = Utilde.F{j}/alphatilde;
                else
                    alphatilde = alpha;
                    Utilde = U;
                end
                global mexcompiled;
                
                % Solver non lineaire
                if (param.bNlSolver)
                    [A,bup,sepmNonLinSys] = UpdateSys(u + alpha*U,solver,b);
                    bu = bup - A*u;
                    UAU = Utilde'*A*U;
                    Ubu = Utilde'*bu;
                    % dNewtonErr = norm(alpha*UAU-Ubu)/norm(Ubu);
                    dNewtonErr = norm(U-Up)/norm(Up);
                    Up=U;
                else
                    if ~mexcompiled
                        for ii=1:A.m
                            AU.F{ii,j} = A.F{ii,j}*U.F{j};
                            UAU.F{ii,j} = Utilde.F{j}'*AU.F{ii,j};
                        end
                        for ii=1:bu.m
                            Ubu.F{ii,j} = Utilde.F{j}'*Ubu.F{ii,j};
                        end
                    else
                        AU.F(:,j) = multiplyF(A.F(:,j),{U.F{j}});
                        UAU.F(:,j) = multiplyF({Utilde.F{j}'},AU.F(:,j));
                        Ubu.F(:,j) = multiplyF({Utilde.F{j}'},bu.F(:,j));
                    end
                    dNewtonErr = -1;
                end
            end %%%% END WHILE NEWTON
        end %%%% END BOUCLE DIMENSION j
        
        erriter(kkk)=abs(alpha-alpha0)/(alpha+alpha0);
        result.alpha{i}(kkk)=alpha;
        result.alphatilde{i}(kkk)=alphatilde;
        
        if param.display && ~(param.node)
            fprintf('%d ',param.node);
            fprintf('iteration #%d - stagnation = %d\n',kkk,erriter(kkk))
        elseif param.display && param.node && kkk==1
            fprintf('%d',param.node);for cpt=1:param.depth, fprintf('      ');end
            fprintf('iteration #%d- %.1d ',kkk,erriter(kkk))
        elseif param.display && param.node && kkk~=1
            fprintf('#%d- %.1d ',kkk,erriter(kkk))
        end
        
        if erriter(kkk)<param.itercrit
            break
        end
        alpha0=alpha;
        alpha0tilde=alphatilde;
        
        U0 = U;
        U0tilde=Utilde;
        
    end  %%%%  END ITERATION ALTERNEE
    
    if param.display && param.node && kkk~=1
        fprintf('\n')
    end
    
    u = u + alpha*U;
    if param.adjoint
        utilde = utilde + alphatilde*Utilde;
    end
    
    
    if param.update>0 && getm(u)>1  && mod(i,param.updatestep)==0
        if param.adjoint
            [u,utilde]=solve_update(A,b,u,utilde,Ametric,param);
        else
            if (param.bNlSolver)
                [u]=solve_update(A,bup,u,[],Ametric,param);
            else
                [u]=solve_update(A,b,u,[],Ametric,param);
            end
        end
        
        if (param.bNlSolver)
            bu = b - sepmNonLinSys*u;
        else
            bu = b - A*u;
        end
    else
        if (param.bNlSolver)
            bu = b - sepmNonLinSys*u;
        else
            bu = b - A*u;
        end
    end
    
    if param.updatetucker && getm(u)>1 && mod(i,param.updatestep)==0 %&& param.update==0
        if param.display; fprintf('%d ',param.node);for cpt=1:param.depth, fprintf('      ');end; fprintf(' updatetucker\n');end
        
        if param.adjoint && param.updateadjoint
            [u,utilde] = solve_updatetucker(A,b,u,utilde,param);
        else
            u = solve_updatetucker(A,b,u,[],param);
        end
        
        bu = b - A*u;
    end
    
    
    if param.alphaupdate && mod(i,param.updatestep)==0
        if param.display;fprintf('%d ',param.node);for cpt=1:param.depth, fprintf('      ');end;fprintf(' alphaupdate\n');end
        if param.adjoint && param.updateadjoint
            u = solve_alphaupdate(A,b,u,utilde);
        else
            u = solve_alphaupdate(A,b,u);
        end
        bu = b - A*u;
    end
    
    
    if param.righthandSD && mod(i,param.righthandSDstep)==0
        mm = getm(bu);
        if getdim(A)>2
            bu = multisvd(bu,'tol',param.tol/5,'maxorder',max(size(A)),'display',false);
            if param.display;fprintf('%d ',param.node);for cpt=1:param.depth, fprintf('      ');end;fprintf('SVD of right-hand side : %d -> %d\n',mm,getm(bu));end
        else
            bu=svd(bu,param.tol/5);
            if param.display;fprintf('%d ',param.node);for cpt=1:param.depth, fprintf('      ');end;fprintf('SVD of right-hand side : %d -> %d\n',mm,getm(bu)); end
        end
    end
    
    if param.storeiter==1
        result.u{i}=u;
    end
    
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
    if param.display
        fprintf('%d',param.node);for cpt=1:param.depth, fprintf('      ');end
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
        fprintf('%d',param.node);for cpt=1:param.depth, fprintf('      ');end
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

function [sepmUpdatedSys,sepvUpdatedeSrc,sepmNonLinSys] =...
    UpdateSys(sepvSol,solvCurrent,sepvSrc)
Params = getparam(solvCurrent);
iRankSvd = 1;
iNbUpSvd = 2;
dErrSvd = 1e-6;
iDim = sepvSol.dim;
bDisp = false;

% Actualisation du systeme, terme de reaction, formulation continue
sepvUWeight = sepvSol.*sepvSol;
solvSvd = setparam(solvCurrent,'maxorder',iRankSvd,'update',iNbUpSvd,...
    'errorindicator','residual','tol',dErrSvd,'display',bDisp);
sepvUWeight = multisvd(sepvUWeight,solvSvd);
cmdvWeight = sepvUWeight.F;
cmblinfD00 = cellfun(@(x)BILINFORM(0,0,x,0),cmdvWeight,'UniformOutput',false);

% Formulation discrete
cmdmD00 = cell(size(cmdvWeight));
for m=1:size(cmdvWeight,1)
    for d=1:iDim
        cmdmD00{m,d} = cmblinfD00{m,d}{Params.modM1d}(:,:);
    end
end

% Assemblage des dimensions
cmdmSys = [Params.cmdmLinOp;cmdmD00];

% Passage aux sepmatrix
dvAlpha = [ones(1,iDim),3*sepvUWeight.alpha]; % derivation de l'operateur
sepmUpdatedSys = SEPMATRIX(cmdmSys,dvAlpha);
dvAlpha = 2*sepvUWeight.alpha; % derivation de l'operateur - 1
sepmPerturbSys = SEPMATRIX(cmdmD00,dvAlpha);
sepvUpdatedeSrc = sepvSrc + sepmPerturbSys*sepvSol;
% sepvUpdatedeSrc = multisvd(sepvUpdatedeSrc,solvSvd);
dvAlpha = [ones(1,iDim),sepvUWeight.alpha]; % operateur non lineaire
sepmNonLinSys = SEPMATRIX(cmdmSys,dvAlpha);

% function M = reduce(A,j)
%
% Z = A.F(:,j);
% A.F(:,j) = [];
% a = A.alpha.*prod(cell2mat(A.F),2)';
% M = a(1)*Z{1};
% for i=2:A.m
% M = M + a(i)*Z{i};
% end
%
% return


%
%
% function M = timesblock(A)
%
% i=1;
% Mi = A.F{i,1};
% for k=2:A.dim
% Mi = Mi.*A.F{i,k};
% end
% M = Mi*A.alpha(i);
%
% for i=2:A.m
% Mi = A.F{i,1};
% for k=2:A.dim
% Mi = Mi.*A.F{i,k};
% end
% M = M + Mi*A.alpha(i);
% end
%
% return






%
% function M = assembleblock(A,j)
%
% Z = A.F(:,j);
% A.F(:,j) = [];
% n = size(Z{1},1);
% nb = size(A.F{1},1);
% B = cell(A.m,1);
% for i=1:A.m
% B{i} = A.F{i,1};
% for k=2:size(A.F,2)
% B{i}=B{i}.*A.F{i,k};
% end
% end
%
%
%
% if size(Z{1},2)>1
%
% nz = zeros(1,length(Z));
% for i=1:length(Z)
%    nz(i)=nnz(Z{i});
% end
% %M = spalloc(n*nb,n*nb,sum(nz)*nb^2);
% %M = zeros(n*nb);
%
% I = zeros(max(nz),A.m,nb,nb);
% J = I;
% V = I;
% for ii=1:nb
% for jj=1:nb
% for i=1:A.m
% [It,Jt,Vt] = find((A.alpha(i)*B{i}(ii,jj))*Z{i});
% I(1:length(It),i,ii,jj)=It + (ii-1)*n;
% J(1:length(It),i,ii,jj)=Jt + (jj-1)*n;
% V(1:length(It),i,ii,jj)=Vt;
% end
%
% %I(1:length(It),ii,jj) = It;
% %J(1:length(It),ii,jj) = Jt;
% %V(1:length(It),ii,jj) = Vt;
% end
% end
% rep = find(V(:));
% I = I(rep);
% J = J(rep);
% V = V(rep);
% M = sparse(I(:),J(:),V(:),n*nb,n*nb);
%
% else
% M = zeros(nb*n,1);
% for ii=1:nb
% I = (ii-1)*n+(1:n);
% MI = (A.alpha(1)*B{1}(ii))*Z{1};
% for i=2:A.m
% MI = MI + (A.alpha(i)*B{i}(ii))*Z{i};
% end
% M(I)=MI;
% end
% end
%
%
%
% return
%
