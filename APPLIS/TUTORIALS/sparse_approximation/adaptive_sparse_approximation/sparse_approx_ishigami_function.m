%% Adaptive sparse polynomial approximation - Ishigami function %%
%%--------------------------------------------------------------%%
% [Ishigami, Homma, 1990, ISUMA]
% [Saltelli, Chan, Scott, 2000, Wiley]
% [Sudret, 2008, RESS]
% [Blatman, Sudret, 2010, PEM]
% [Blatman, Sudret, 2010, RESS]
% [Blatman, Sudret, 2011, JCP]

% clc
clearvars
close all
rng('default');

%% Random variables
M = 3; % number of random variables
rv = RVUNIFORM(-pi,pi);
RV = RANDVARS(repmat({rv},1,M));

%% Ishigami function
% y = sin(x_1) + a*sin(x_2)^2 + b*(x_3)^4*sin(x_1)
a = 7;
b = 0.1;
fun = @(x) (sin(x(:,1)) + a.*sin(x(:,2)).^2 + b.*x(:,3).^4.*sin(x(:,1)))';

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
opts.tol = 1e-4; % prescribed tolerance for cross-validation error
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
fprintf('error = %e (cross-validation error)\n',err)
fprintf('elapsed time = %f s\n',time)

%% Display evolution of multi-index set
dim = 1:3;
video_indices(PC_seq,'dim',dim)

%% Display evolution of cross-validation error indicator, dimension of stochastic space and number of samples w.r.t. number of iterations
% plot_adaptive_algorithm(err_seq,PC_seq,N_seq);

%% Display evolution of cross-validation error indicator w.r.t. number of samples
% plot_cv_error_indicator_vs_nb_samples(err_seq,N_seq,'legend',false);

%% Display evolution of cross-validation error indicator w.r.t. dimension of stochastic space
% plot_cv_error_indicator_vs_dim_stochastic_space(err_seq,PC_seq,'legend',false);

%% Display multi-index set
dim = 1:3;
plot_multi_index_set(PC,'dim',dim,'legend',false)

%% Quantities of interest : mean, variance, Sobol indices
% Analytical exact values
anal.mean = a/2;
anal.var = 1/2 + a^2/8 + b*pi^4/5 + b^2*pi^8/18;
anal.S1 = (1/2*(1 + b*pi^4/5)^2)/anal.var;
anal.S2 = (a^2/8)/anal.var;
anal.S3 = 0;
anal.S12 = 0;
anal.S13 = (b^2*pi^8*(8/225))/anal.var;
anal.S23 = 0;
anal.S123 = 0;
anal.S1T = anal.S1 + anal.S12 + anal.S13 + anal.S123;
anal.S2T = anal.S2 + anal.S12 + anal.S23 + anal.S123;
anal.S3T = anal.S3 + anal.S13 + anal.S23 + anal.S123;
% Numerical approximate values
num.mean = mean(u);
num.var = variance(u);
num.S1 = sobol_indices(u,1);
num.S2 = sobol_indices(u,2);
num.S3 = sobol_indices(u,3);
num.S12 = sobol_indices_group(u,[1,2]) - num.S1 - num.S2;
num.S13 = sobol_indices_group(u,[1,3]) - num.S1 - num.S3;
num.S23 = sobol_indices_group(u,[2,3]) - num.S2 - num.S3;
num.S123 = sobol_indices_group(u,[1,2,3]) - num.S1 - num.S2 - num.S3 - num.S12 - num.S13 - num.S23;
num.S1T = num.S1 + num.S12 + num.S13 + num.S123;
num.S2T = num.S2 + num.S12 + num.S23 + num.S123;
num.S3T = num.S3 + num.S13 + num.S23 + num.S123;
% Comparative table
fanal = '%10.5f';
fnum = '%9.5f';
ferr = '%14.4e';
fprintf('\n')
disp('+-------------------+------------+-----------+----------------+')
disp('| Quantity \ Value  | Analytical | Numerical | Relative error |')
disp('+-------------------+------------+-----------+----------------+')
fprintf(['| Mean E            | ' fanal ' | ' fnum ' | ' ferr ' |\n'],anal.mean,num.mean,abs((anal.mean-num.mean)/anal.mean))
fprintf(['| Variance V        | ' fanal ' | ' fnum ' | ' ferr ' |\n'],anal.var,num.var,abs((anal.var-num.var)/anal.var))
fprintf(['| Sobol index S_1   | ' fanal ' | ' fnum ' | ' ferr ' |\n'],anal.S1,num.S1,abs((anal.S1-num.S1)/anal.S1))
fprintf(['| Sobol index S_2   | ' fanal ' | ' fnum ' | ' ferr ' |\n'],anal.S2,num.S2,abs((anal.S2-num.S2)/anal.S2))
fprintf(['| Sobol index S_3   | ' fanal ' | ' fnum ' |                |\n'],anal.S3,num.S3)
fprintf(['| Sobol index S_12  | ' fanal ' | ' fnum ' |                |\n'],anal.S12,num.S12)
fprintf(['| Sobol index S_13  | ' fanal ' | ' fnum ' | ' ferr ' |\n'],anal.S13,num.S13,abs((anal.S13-num.S13)/anal.S13))
fprintf(['| Sobol index S_23  | ' fanal ' | ' fnum ' |                |\n'],anal.S23,num.S23)
fprintf(['| Sobol index S_123 | ' fanal ' | ' fnum ' |                |\n'],anal.S123,num.S123)
fprintf(['| Sobol index S_1^T | ' fanal ' | ' fnum ' | ' ferr ' |\n'],anal.S1T,num.S1T,abs((anal.S1T-num.S1T)/anal.S1T))
fprintf(['| Sobol index S_2^T | ' fanal ' | ' fnum ' | ' ferr ' |\n'],anal.S2T,num.S2T,abs((anal.S2T-num.S2T)/anal.S2T))
fprintf(['| Sobol index S_3^T | ' fanal ' | ' fnum ' | ' ferr ' |\n'],anal.S3T,num.S3T,abs((anal.S3T-num.S3T)/anal.S3T))
disp('+-------------------+------------+-----------+----------------+')
