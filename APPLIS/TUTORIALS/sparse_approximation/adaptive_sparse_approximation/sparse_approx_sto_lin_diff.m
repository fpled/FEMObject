%% Monoscale stochastic linear diffusion %%
%%---------------------------------------%%

% clc
clearvars
close all

% Parallel computing
myparallel('start');

%% Domains and meshes

D = DOMAIN(2,[0.0,0.0],[1.0,1.0]);

nbelem = [20,20];
system.S = build_model(D,'nbelem',nbelem);
% cl = 0.05;
% system.S = build_model(D,'cl',cl,'filename','gmsh_domain');

%% Random variables

rv = RVUNIFORM(0,1);
RV = RANDVARS(rv);
[X,PC] = PCMODEL(RV,'order',1,'pcg','typebase',1);

%% Materials

% Linear diffusion coefficient
% K(xi) = 1 + xi
K = ones(1,1,PC) + X{1};
mat = FOUR_ISOT('k',K); % uniform value
system.S = setmaterial(system.S,mat);

%% Dirichlet boundary conditions

system.S = final(system.S);
system.S = addcl(system.S,[]);

%% Stiffness matrices and sollicitation vectors

% Source term
f = 100;

if israndom(system.S)
    system.A = [];
else
    system.A = calc_rigi(system.S);
end

if israndom(f)
    system.b = [];
else
    system.b = bodyload(system.S,[],'QN',f);
end

%% Resolution

initPC = POLYCHAOS(RV,0,'typebase',1);

method = METHOD('type','leastsquares','display',true,'displayiter',true,...
    'basis','adaptive','initPC',initPC,'maxcoeff',Inf,...
    'algorithm','RMS','bulkparam',0.5,...
    'sampling','adaptive','initsample',2,'addsample',0.1,'maxsample',Inf,...
    'regul','','cv','leaveout','k',10,...
    'tol',1e-12,'tolstagn',1e-1,'toloverfit',1.1,'correction',false,...
    'decompKL',false,'tolKL',1e-12,'cvKL','leaveout','kKL',10);

fun = @(xi) solve_system(calc_system(randomeval_system(system,xi)));
[u,result] = solve_random(method,fun);
PC = getPC(u);

%% Display domains and meshes

plot_domain(D);
% plot_partition(system.S,'legend',false);
plot_model(system.S,'legend',false);

%% Display multi-index set

plot_multi_index_set(PC,'legend',false)

%% Display evolution of multi-index set

if isfield(result,'PC_seq')
    video_indices(result.PC_seq)
end

%% Display evolution of cross-validation error indicator, dimension of stochastic space and number of samples w.r.t. number of iterations

if isfield(result,{'cv_error_indicator_seq','PC_seq','N_seq'})
    plot_adaptive_algorithm(result.cv_error_indicator_seq,result.PC_seq,result.N_seq);
end

%% Display statistical outputs of solution

% plot_stats(system.S,u);

plot_mean(system.S,u);
plot_var(system.S,u);
plot_std(system.S,u);

M = getM(PC);
for m=1:M
    plot_sobol_indices(system.S,u,m);
    plot_sensitivity_indices_max_var(system.S,u,m);
end

%% Display random evaluations of solution

% nbSamples = 3;
% for s=1:nbSamples
%     xi = random(RANDVARS(PC));
%     u_xi = randomeval(u,xi);
%     S_xi = randomeval(system.S,xi);
%     plotSolution(S_xi,u_xi);
% end

myparallel('stop');
