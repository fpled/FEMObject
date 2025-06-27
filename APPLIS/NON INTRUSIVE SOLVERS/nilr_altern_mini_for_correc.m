function [dV,dLambda,result] = nilr_altern_mini_for_correc(pw_J,pw_R,sz_v,V,Lambda,opts)
% function [dV,dLambda,result] = nilr_altern_mini_for_correc(pw_J,pw_R,sz_v,,V,Lambda,opts)
%
% Compute a low-rank correction dV*dLambda'of the approximation of the
% solution u satisfying  
% pw_J(p,u(p)) = min pw_J(p,v), for all p.
% 
% The low-rank approximation of u is V*Lambda' in an algebraic context, the
% new approximation at the end should be [V dV]*[Lambda dLambda]'
%
% pw_R is the residual of the equation, that is 
% pw_R(p,v) = -\nabla_v pw_J(p,v)
%
%
% **** Options concerning the solvers are
% rank : rank of the low-rank approximation (default 10)
% maxiter : number of iterations of the alternating minimization algorithm 
%           (default 10)
% tol : stagnation criterion used for the AMA (default 1e-8)
% fminuncopt : options for the fminunc function (default is
% optimset('GradObj','on','LargeScale','off','MaxIter',200,'TolFun',100*eps,'TolX',100*eps))
% 
%
% **** Options concerning the integration are
% PC : automatically choose quadrature points (default quad_order = 1);
% PC + quad_order : change the quad_order
% PC + gp : automatically evaluates the polynomial at Gauss points
% gp + polyval : manually supply the Gauss points and the evaluations of
% polynomial at Gauss points

timetemp = tic;
opts = defaultopts(opts);

sz_lam = size(opts.polyval,2);

if isfield(opts,'lambda_init') && ~isempty(opts.lambda_init)
    dLambda = opts.lambda_init;
else
    dLambda = rand(sz_lam,opts.rank);
end

if isfield(opts,'v_init') && ~isempty(opts.v_init)
    dV = opts.v_init;
else
    dV = rand(sz_v,opts.rank);
end

result.error = zeros(1,opts.maxiter);

gp = opts.gp;
polyval = opts.polyval;
r = opts.rank;

i = 1;
erriter = 1;
while (i <= opts.maxiter) && (erriter > opts.tol)
    fprintf('AMA : Iter %d \n',i)
    dLambda0 = dLambda;
    dV0 = dV;
    if r > 1
        [dV,R] = qr(dV,0);
        dLambda = dLambda*R';
    else
        nv = norm(dV);
        dV = dV/nv;
        dLambda = nv*dLambda;
    end
    %% Minimization on Lambda
    fprintf('  mini Lambda \n')
    JLam = @(lam) nilr_int_cost_fun([V dV],[Lambda lam],pw_J,gp,polyval);
    GLam = @(lam) -eval_RLambda(V,dV,Lambda,lam,pw_R,gp,polyval);
    dLambda = fminunc(@(lam) myfun(lam,JLam,GLam,sz_lam,r),dLambda(:),opts.fminuncopt);
    dLambda = reshape(dLambda,sz_lam,r);
    if r > 1;
        [dLambda,R] = qr(dLambda,0);
        dV = dV*R';
    else
        nlambda = norm(dLambda);
        dLambda = dLambda/nlambda;
        dV = nlambda*dV;
    end
    
    %% Minimization on V
    fprintf('  mini V\n')
    JV = @(v) nilr_int_cost_fun([V v],[Lambda dLambda],pw_J,gp,polyval);
    GV = @(v) -eval_RV(V,v,Lambda,dLambda,pw_R,gp,polyval);
    dV = fminunc(@(v) myfun(v,JV,GV,sz_v,r),dV(:),opts.fminuncopt);
    dV = reshape(dV,sz_v,opts.rank);
    
    %%
    erriter = norm((dV*dLambda')-(dV0*dLambda0'),'fro')/...
        norm(dV0*dLambda0','fro');
    result.error(i) = erriter;
    fprintf('\t\t err %d\n',erriter);
    i = i+1;
end
result.time = toc(timetemp);

end

function [RLam] = eval_RLambda(V,dV,Lambda,dLambda,pw_R,gp,polyval)
Rz = [V dV]*(polyval*[Lambda dLambda])';
for z = 1:gp.nbgauss
    Rz(:,z) = pw_R(gp.coord(z,:),Rz(:,z));
end

PsiR = gp.w(1)*polyval(1,:)'*Rz(:,1)';
for z = 2:gp.nbgauss
    PsiR = PsiR + gp.w(z)*polyval(z,:)'*Rz(:,z)';
end

RLam = PsiR*dV;
RLam = RLam(:);
end

function [RV] = eval_RV(V,dV,Lambda,dLambda,pw_R,gp,polyval)

Rz = [V dV]*(polyval*[Lambda dLambda])';
for z = 1:gp.nbgauss
    Rz(:,z) = pw_R(gp.coord(z,:),Rz(:,z));
end

RPsi = gp.w(1)*Rz(:,1)*polyval(1,:);
for z = 2:gp.nbgauss
    RPsi = RPsi + gp.w(z)*Rz(:,z)*polyval(z,:);
end

RV = RPsi*dLambda;
RV = RV(:);
end


function [J,G] = myfun(x,J,G,n,r)
x = reshape(x,n,r);
J = J(x);
G = G(x);
G = G(:);
end

function opts = defaultopts(opts)
    % Parameters for alternating minim
    if ~isfield(opts,'rank');opts.rank = 10; end
    if ~isfield(opts,'maxiter');opts.maxiter = 10; end
    if ~isfield(opts,'tol');opts.tol = 1e-8; end
    
    % Parameters for local optimization
    if ~isfield(opts,'fminuncopt')
        opts.fminuncopt = optimset('GradObj','on','LargeScale','off',...
            'MaxIter',200,'TolFun',100*eps,'TolX',100*eps);...
    end
    
    % Parameters for integration
    if isfield(opts,'PC')
        if ~isfield(opts,'gp')    
            if isfield(opts,'quad_order')
                opts.gp = calc_gausspoints(opts.PC,opts.quad_order);
            else
                opts.gp = calc_gausspoints(opts.PC,getorder(opts.PC)+1);
            end
        end
        if ~isfield(opts,'polyval')
            opts.polyval = full(polyval(opts.PC,opts.gp.coord));
        end
    else
        if ~isfield(opts,'gp')
            error('Gauss points are missing')
        end
        if ~isfield(opts,'polyval')
            error('Evaluations of polynomial at Gauss points are missing')
        end
    end
    
    
    if opts.rank > size(opts.polyval,2)
        error('opts.rank must be inferior or equal to %d',size(opts.polyval,2))
    end
end


