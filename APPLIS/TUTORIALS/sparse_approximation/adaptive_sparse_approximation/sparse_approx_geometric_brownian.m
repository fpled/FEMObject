%% Adaptive sparse polynomial approximation - Geometric brownian %%
%%---------------------------------------------------------------%%

% clc
clearvars
close all
rng('default');

%% Random variables
M = 10; % number of random variables
rv = RVNORMAL(0,1);
RV = RANDVARS(repmat({rv},1,M));

%% Geometric brownian
c = -1;
s = 0.5;
T = 1;
m = 101; % output size
n = m-1;
fun = @(x) geometric_brownian(x,c,s,T,n)';

%% Resolution using adaptive sparse approximation and least-squares minimization

% Polynomial chaos basis
p = 0; % (initial) order of PC expansion
PC = POLYCHAOS(RV,p,'typebase',1); % (initial) PC basis
opts = struct();
opts.basis = 'adaptive'; % construction of PC basis ('fixed' or 'adaptive')
opts.maxcoeff = Inf; % maximal number of unknown PC expansion coefficients
opts.algorithm = 'RMS'; % adaptive algorithm for the construction of a nested sequence of finite monotone/lower multi-index sets
% 'TP' or 'PD':  isotropic Tensor Product (or Partial Degree) polynomial space
%        multidimensional space of polynomials of partial degree less or equal to p (in each variable)
%        update the partial degree by 1 in each dimension at each iteration
% 'TD':  isotropic Total Degree polynomial space
%        multidimensional space of polynomials of total degree less or equal to p
%        update the total degree by 1 at each iteration
% 'MS':  Margin Search strategy
%        add a smallest monotone/lower subset S_n of the margin M_n of a given monotone/lower set A_n
%        for which energy(S_n)>=bulkparam*energy(M_n), where bulkparam is a bulk parameter
% 'RMS': Reduced Margin Search strategy
%        add a smallest monotone/lower subset S_n of the reduced margin M_n of a given monotone/lower set A_n
%        for which energy(S_n)>=bulkparam*energy(M_n), where bulkparam is a bulk parameter
opts.bulkparam = 0.5; % bulk parameter in (0,1) such that energy(S_n)>=bulkparam*energy(M_n)
% bulkparam = 1 selects all multi-indices in the (reduced) margin M_n
% bulkparam = 0 selects the multi-index in the (reduced) margin M_n
%               corresponding to the maximal norm of the expansion coefficients
funtr = @(x) fun(transfer(RANDVARS(PC),RANDVARS(RV),x));

% Sampling
N = 1; % (initial) number of samples (regression points)
opts.sampling = 'adaptive'; % sampling strategy ('fixed' or 'adaptive')
opts.addsample = 0.1; % percentage of additional samples if 0 < addsample < 1
                      % number of additional samples if addsample > 1
opts.maxsample = Inf; % maximal number of samples

% Regularization
regul = ''; % type of regularization ('' or 'l0' or 'l1')

% Cross-validation
cv = 'leaveout'; % type of cross-validation procedure ('leaveout' or 'kfold')
k = 10; % number of folds (only for k-fold cross-validation procedure)
opts.tol = 1e-3; % prescribed tolerance for cross-validation error
opts.tolstagn = 5e-2; % prescribed stagnation tolerance for cross-validation error
opts.toloverfit = 1.1; % prescribed tolerance to detect overfitting for cross-validation error such that err>=toloverfit*err_old
opts.correction = true; % correction for cross-validation error

% Least-squares minimization
t = tic;
[u,err,x,PC_seq,err_seq,x_seq] = decompmatrix_leastsquares(PC,N,funtr,1,regul,[],cv,k,opts,'displayiter');
N_seq = cellfun(@(x) size(x,1),x_seq);
PC = getPC(u);
time = toc(t);

%% Results
disp(' ')
disp(['M = ' num2str(getM(PC)) ' random variables'])
disp(['p = ' num2str(getorder(PC)) ' (order of PC expansion)'])
disp(['P = ' num2str(length(PC)) ' unknown PC expansion coefficients'])
disp(['N = ' num2str(size(x,1)) ' samples'])
disp(['I = ' num2str(size(getindices(PC),1)) ' multi-indices']);
% disp('Set of multi-indices = '); % P-by-(M+1) matrix
% disp(num2str(getindices(PC)));
disp(['eta = ' num2str(getSparsityRatio(u)) ' (sparsity index or ratio)'])
fprintf('error = %e (cross-validation error)\n',norm(err))
fprintf('elapsed time = %f s\n',time)

%% Display random evaluations of Brownian motion
xi = cell2mat(random(RV));
Xref = fun(xi);
Xpc = randomeval(u,xi);
plot_geometric_brownian(Xref,Xpc);

%% Display evolution of multi-index set
video_indices(PC_seq,'dim',[1 3 5])
video_indices(PC_seq,'dim',[1 7 10])

%% Display evolution of cross-validation error indicator, dimension of stochastic space and number of samples w.r.t. number of iterations
% plot_adaptive_algorithm(err_seq,PC_seq,N_seq);

%% Display evolution of cross-validation error indicator w.r.t. number of samples
% plot_cv_error_indicator_vs_nb_samples(err_seq,N_seq);

%% Display evolution of cross-validation error indicator w.r.t. dimension of stochastic space
% plot_cv_error_indicator_vs_dim_stochastic_space(err_seq,PC_seq);

%% Display multi-index set
dim = [1 3 5];
plot_multi_index_set(PC,'dim',dim,'legend',false)

dim = [1 7 10];
plot_multi_index_set(PC,'dim',dim,'legend',false)
