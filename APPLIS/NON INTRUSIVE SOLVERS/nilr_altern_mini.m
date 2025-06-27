function [V,Lambda,result] = nilr_altern_mini(pw_J,pw_R,sz_v,opts)
% function [V,Lambda,result] = nilr_altern_mini(pw_J,pw_R,sz_v,opts)
%
% Compute a low-rank approximation V*Lambda'of the solution u satisfying 
% pw_J(p,u(p)) = min pw_J(p,v), for all p.
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
    Lambda = opts.lambda_init;
else
    Lambda = rand(sz_lam,opts.rank);
end

if isfield(opts,'v_init') && ~isempty(opts.v_init)
    V = opts.v_init;
else
    V = rand(sz_v,opts.rank);
end

result.error = zeros(1,opts.maxiter);

gp = opts.gp;
polyval = opts.polyval;
r = opts.rank;

i = 1;
erriter = 1;
while (i <= opts.maxiter) && (erriter > opts.tol)
    fprintf('AMA : Iter %d \n',i)
    Lambda0 = Lambda;
    V0 = V;
    if r > 1
        [V,R] = qr(V,0);
        Lambda = Lambda*R';
    else
        nv = norm(V);
        V = V/nv;
        Lambda = nv*Lambda;
    end
    %% Minimization on Lambda
    fprintf('  mini Lambda \n')
    JLam = @(lam) nilr_int_cost_fun(V,lam,pw_J,gp,polyval);
    GLam = @(lam) -eval_RLambda(V,lam,pw_R,gp,polyval);
    Lambda = fminunc(@(lam) myfun(lam,JLam,GLam,sz_lam,r),Lambda(:),opts.fminuncopt);
    Lambda = reshape(Lambda,sz_lam,r);
    if r > 1;
        [Lambda,R] = qr(Lambda,0);
        V = V*R';
    else
        nlambda = norm(Lambda);
        Lambda = Lambda/nlambda;
        V = nlambda*V;
    end
    
    %% Minimization on V
    fprintf('  mini V\n')
    JV = @(v) nilr_int_cost_fun(v,Lambda,pw_J,gp,polyval);
    GV = @(v) -eval_RV(v,Lambda,pw_R,gp,polyval);
    V = fminunc(@(v) myfun(v,JV,GV,sz_v,r),V(:),opts.fminuncopt);
    V = reshape(V,sz_v,opts.rank);
    
    %%
    erriter = norm((V*Lambda')-(V0*Lambda0'),'fro')/...
        norm(V0*Lambda0','fro');
    result.error(i) = erriter;
    fprintf('\t\t err %d\n',erriter);
    i = i+1;
end
result.time = toc(timetemp);

end

function [RLam] = eval_RLambda(V,Lambda,pw_R,gp,polyval)

Rz = V*(polyval*Lambda)';
for z = 1:gp.nbgauss
    Rz(:,z) = pw_R(gp.coord(z,:),Rz(:,z));
end

PsiR = gp.w(1)*polyval(1,:)'*Rz(:,1)';
for z = 2:gp.nbgauss
    PsiR = PsiR + gp.w(z)*polyval(z,:)'*Rz(:,z)';
end

RLam = PsiR*V;
RLam = RLam(:);

end

function [RV] = eval_RV(V,Lambda,pw_R,gp,polyval)

Rz = V*(polyval*Lambda)';
for z = 1:gp.nbgauss
    Rz(:,z) = pw_R(gp.coord(z,:),Rz(:,z));
end

RPsi = gp.w(1)*Rz(:,1)*polyval(1,:);
for z = 2:gp.nbgauss
    RPsi = RPsi + gp.w(z)*Rz(:,z)*polyval(z,:);
end

RV = RPsi*Lambda;
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


