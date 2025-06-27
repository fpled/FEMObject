clear all

n = 40;
Sx = mesh(LIGNE(0,1),n);
x = getcoord(getnode(Sx));
Sx = createddlnode(Sx,DDL('U'));
Sx = addcl(Sx,POINT([0;1]),'U');
ax = BILINFORM(1,1);
fx = LINFORM(0,-1);


Ax = ax{Sx}(:,:);

Fx = fx{Sx}(:);

%% Obstacle
xf = freevector(Sx,x);
gtx = @(t,x) t*max(sin(3*pi*x),0) + (t-1)*min(sin(3*pi*x),0);

gt = @(t) gtx(t,x);  
gtf = @(t) gtx(t,xf);

%% Point-wise parametric maps
epsilon = 1000;

fun_param = @(t,x) -max(gtf(t)-x,0);
funtang_param = @(t,x) double(x<double(gtf(t)));

J_param = @(t,w) 0.5*w'*Ax*w-w'*Fx;
Jpen_param = @(t,w) J_param(t,w) + 0.5*epsilon*(fun_param(t,w)'*fun_param(t,w));
calcAu_param = @(t,w) Ax*w+epsilon*fun_param(t,w);
calcAtang_param = @(t,w) Ax+diag(epsilon*funtang_param(t,w));

pw_residual = @(t,w) Fx-calcAu_param(t,w);


%% L2 projection of the solution
RV = RVUNIFORM(0,1);

ord = 1;
[X,PC] = PCMODEL(RV,'order',ord,'fedim',1,'femesh',{20});

numgp = ord+1;
tol = 1e-12;
maxiter = 200;

NS = NEWTONSOLVER('type','full','display',false,'tol',tol,'maxiter',maxiter);
eval_u = @(t) solve(NS,Fx,@(w) calcAu_param(t,w),@(w) calcAtang_param(t,w));
c=clock;
upc = decompmatrixiter(PC,numgp,[],eval_u);
etime(clock,c)

%%
D = svd(double(upc));

err_svd = sqrt(1-cumsum(D.^2)/sum(D.^2));

figure(555)
clf
semilogy(err_svd)

%%
T = 0:0.02:1;
upctot = double(randomeval(upc,T'));

figure(1)
clf
surf(upctot)

%% 
N_samples = 2000;
[upc_samples,xi_samples] = random(upc,N_samples);
upc_samples = double(upc_samples);

uref_samples = zeros(n-1,N_samples);

for i = 1:N_samples
    uref_samples(:,i) = eval_u(xi_samples(i)); 
end
%%
err = zeros(N_samples,1);
for i = 1:N_samples
    err(i) = norm(uref_samples(:,i)-upc_samples(:,i))^2;
end
RMSE = sqrt(mean(err))

%% TEST LOW RANK
sz_v = getnbddlfree(Sx);

clear opts
opts.rank = 5;
opts.maxiter = 10;
opts.tol = 1e-8;
opts.fminuncopt = optimset('GradObj','on','LargeScale','off',...
            'MaxIter',200,'TolFun',100*eps,'TolX',100*eps);
opts.PC = PC;
[V,Lambda,result] = nilr_altern_mini(Jpen_param,pw_residual,sz_v,opts);
result
%% DIRECT MINIMIZATION

maxrank = 9;
sz_lambda = getn(PC);

clear opts
opts.rank = 5;
opts.maxiter = 10;
opts.tol = 1e-8;
opts.fminuncopt = optimset('GradObj','on','LargeScale','off',...
            'MaxIter',200,'TolFun',100*eps,'TolX',100*eps);
opts.PC = PC;

if isfield(opts,'v_init');opts = rmfield(opts,'v_init');end
if isfield(opts,'lambda_init');opts = rmfield(opts,'lambda_init');end

V_ama = cell(maxrank,1);
Lambda_ama = V_ama;
result_ama = V_ama;

time_ama = zeros(maxrank,1);

for r = 1:maxrank
    timetemp = tic;
    opts.rank = r;
    if r > 1
        opts.v_init = zeros(sz_v,r);
        opts.v_init(:,1:r-1) = V_ama{r-1};
        opts.lambda_init = zeros(sz_lambda,r);
        opts.lambda_init(:,1:r-1) = Lambda_ama{r-1};
    end
    [V_ama{r},Lambda_ama{r},result_ama{r}] = nilr_altern_mini(Jpen_param,pw_residual,sz_v,opts);
    time_ama(r) = toc(timetemp);
end

%%
RMSE_ama = zeros(maxrank,1);

for r = 1:maxrank
    u_ama = PCMATRIX(V_ama{r}*Lambda_ama{r}',[sz_v 1],PC);
    u_ama_samples = double(randomeval(u_ama,xi_samples));
    err = zeros(N_samples,1);
    for i = 1:N_samples
        err(i) = norm(uref_samples(:,i)-u_ama_samples(:,i))^2;
    end
    RMSE_ama(r) = sqrt(mean(err));
end

RMSE_ama

%% GREEDY APPROXIMATION
sz_lambda = getn(PC);

clear opts
opts.rank = 1;
opts.maxiter = 10;
opts.tol = 1e-8;
opts.fminuncopt = optimset('GradObj','on','LargeScale','off',...
            'MaxIter',200,'TolFun',100*eps,'TolX',100*eps);
opts.PC = PC;
if isfield(opts,'v_init');opts = rmfield(opts,'v_init');end
if isfield(opts,'lambda_init');opts = rmfield(opts,'lambda_init');end

opts.maxcorrec = 20;
[V_greedy,Lambda_greedy,result_greedy] = nilr_greedy_approx(Jpen_param,pw_residual,sz_v,opts);

%%
RMSE_greedy = zeros(opts.maxcorrec,1);
for i = 1:opts.maxcorrec
    u_greedy = PCMATRIX(V_greedy(:,1:i)*Lambda_greedy(:,1:i)',[sz_v 1],PC);
    u_greedy_samples = double(randomeval(u_greedy,xi_samples));
    err = zeros(N_samples,1);
    for j = 1:N_samples
        err(j) = norm(u_greedy_samples(:,j)-uref_samples(:,j))^2;
    end
    RMSE_greedy(i) = sqrt(mean(err));
end
RMSE_greedy

%%
figure(999)
clf
semilogy(cumsum(time_ama),RMSE_ama,'+-');
hold on
semilogy(cumsum(result_greedy.time),RMSE_greedy,'o-r')
xlabel('Time (s)')
ylabel('RMSE')
legend('Direct min.','Greedy rank 1')

figure(1000)
clf
semilogy(1:maxrank,RMSE_ama,'+-')
hold on
semilogy(1:opts.maxcorrec,RMSE_greedy,'o-r')
xlabel('Rank')
ylabel('RMSE')
legend('Direct min.','Greedy rank 1');

