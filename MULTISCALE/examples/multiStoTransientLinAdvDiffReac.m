%% Multiscale stochastic transient linear advection-diffusion-reaction problem %%
%%-----------------------------------------------------------------------------%%

% clc
clearvars
close all
% rng('default');
myparallel('start');

%% Input data
setProblem = true;
directSolver = true;
iterativeSolver = true;
displaySolution = true;

n = 1; % number of patches
filename = ['transientLinAdvDiffReac' num2str(n) 'Patches'];
% for rho = [0.2:0.2:1.8 0.9 1.1]
% close all
% filename = ['transientLinAdvDiffReac' num2str(n) 'PatchesRho' num2str(rho)];
% for tol = 1:3
% close all
% filename = ['transientLinAdvDiffReac' num2str(n) 'PatchesRhoAitkenTol'  num2str(tol)];
pathname = fullfile(getfemobjectoptions('path'),'MULTISCALE',...
    'examples','multiscaleSto',filename);
if ~exist(pathname,'dir')
    mkdir(pathname);
end
fontsize = 16;
formats = {'fig','epsc'};
renderer = 'OpenGL';

%% Problem
if setProblem
    %% Domains and meshes
    % Global
    glob = Global();
    globOut = GlobalOutside();
    
    % Patches
    patches = Patches(n);
    
    D_patch = cell(1,n);
    if n>=1
        D_patch{1} = DOMAIN(2,[0.85,0.40],[1.05,0.80]);
    end
    if n>=2
        D_patch{2} = DOMAIN(2,[0.45,0.20],[0.65,0.60]);
    end
    if n>=3
        D_patch{3} = DOMAIN(2,[0.05,0.40],[0.25,0.80]);
    end
    
    D_patch_star = cell(1,n);
    if n>=1
        D_patch_star{1} = DOMAIN(2,[0.85,0.40],[1.05,0.60]);
    end
    if n>=2
        D_patch_star{2} = DOMAIN(2,[0.45,0.40],[0.65,0.60]);
    end
    if n>=3
        D_patch_star{3} = DOMAIN(2,[0.05,0.40],[0.25,0.60]);
    end
    
    cl1 = 0.02;
    cl2 = 0.04;
    cl0 = 0.02;
    cltip = 0.01;
    clI = 0.02;
    glob.S = gmshCanisterMulti(D_patch,cl1,cl2,cl0,cltip,clI,fullfile(pathname,'gmsh_canister_multi'));
    % glob.S = gmsh2femobject(2,fullfile(pathname,'gmsh_canister_multi.msh'),2);
    
    % nbelem_patch = [10,20];
    % for k=1:n
    %     patches.patches{k}.S = build_model(D_patch{k},'nbelem',nbelem_patch);
    % end
    cl_patch = 0.01;
    for k=1:n
        patches.patches{k}.S = build_model(D_patch{k},'cl',cl_patch,'filename',fullfile(pathname,['gmsh_patch_' num2str(k)]));
        % patches.patches{k}.S = gmsh2femobject(2,fullfile(pathname,['gmsh_patch_' num2str(k) '.msh']),2);
    end
    
    % Partition of global mesh
    glob = partition(glob,D_patch);
    
    %% Random variables
    d = 3*n; % parametric dimension
    r = UniformRandomVariable(0,1);
    rv = RandomVector(r,d);
    
    %% Materials
    % Linear diffusion coefficient
    K_out = 0.01;
    K_patch = cell(1,n);
    K_in = cell(1,n);
    % Thermal capacity
    c_out = 1;
    c_patch = cell(1,n);
    c_in = cell(1,n);
    % Advection velocity
    Sadv = glob.S;
    mat = FOUR_ISOT('k',1);
    mat = setnumber(mat,1);
    Sadv = setmaterial(Sadv,mat);
    P = @(i) POINT(getnode(getridge(Sadv,i)));
    L1 = LIGNE(P(5),P(6));
    L2 = LIGNE(P(15),P(16));
    Sadv = final(Sadv);
    Sadv = addcl(Sadv,P(1),'T',0);
    A = calc_rigi(Sadv);
    b1 = surfload(Sadv,L1,'QN',-1);
    b2 = surfload(Sadv,L2,'QN',1);
    b = b1+b2;
    phi = A\b;
    v = 2*FENODEFIELD(calc_sigma(Sadv,phi,'node'));
    V = getvalue(v);
    V = {{FENODEFIELD(V(:,1)),FENODEFIELD(V(:,2))}};
    V_out = V;
    V_in = cell(1,n);
    for k=1:n
        V_in{k} = V;
    end
    Sadv_patch = cell(1,n);
    for k=1:n
        Sadv_patch{k} = patches.patches{k}.S;
        Sadv_patch{k} = setmaterial(Sadv_patch{k},mat);
        Sadv_patch{k} = final(Sadv_patch{k});
    end
    phi_patch = cell(1,n);
    v_patch = cell(1,n);
    V_patch = cell(1,n);
    for k=1:n
        phi_patch{k} = calcProjection(Sadv_patch{k},Sadv)'*phi;
        v_patch{k} = 2*FENODEFIELD(calc_sigma(Sadv_patch{k},phi_patch{k},'node'));
        V_patch{k} = getvalue(v_patch{k});
        V_patch{k} = {{FENODEFIELD(V_patch{k}(:,1)),FENODEFIELD(V_patch{k}(:,2))}};
    end
    % Linear reaction parameter
    R1_out = 0.1;
    R2_out = 10;
    R_patch = cell(1,n);
    R_in = cell(1,n);
    
    % IntegrationRule
    p = 1;
    bases = cellfun(@(x) PolynomialFunctionalBasis(x,0:p),orthonormalPolynomials(rv),'UniformOutput',false);
    bases = FunctionalBases(bases);
    H = FullTensorProductFunctionalBasis(bases);
    I = gaussIntegrationRule(r,2);
    I = I.tensorize(d);
    
    for k=1:n
        patch = patches.patches{k};
        % K_patch(x,xi) = K_out * (1 + 0.25 * (2 * xi - 1) * f(x))
        % K_in(x)       = K_out
        % c_patch(x,xi) = c_out * (1 + 0.1 * (2 * xi - 1) * f(x))
        % c_in(x)       = c_out
        % R_patch(x,xi) = R1_out * (1 + 0.25 * (2 * xi - 1) * f(x))
        % R_in(x)       = R1_out
        % with f(x) = 1 if ||x-c||_Inf < L
        %           = 0 if ||x-c||_Inf >= L
        L = norm(getsize(D_patch_star{k}))/4;
        c = getcenter(D_patch_star{k});
        f = @(x) distance(x,c,Inf)<L;
        
        fun = @(xi) K_out * (ones(size(xi,1),patch.S.nbnode) + 0.25 * (2 * xi(:,3*k-2) - 1) * double(squeeze(f(patch.S.node)))');
        fun = UserDefinedFunction(fun,d,patch.S.nbnode);
        fun.evaluationAtMultiplePoints = true;
        
        K_patch{k} = FENODEFIELD(H.projection(fun,I));
        
        fun = @(xi) c_out * (ones(size(xi,1),patch.S.nbnode) + 0.1 * (2 * xi(:,3*k-1) - 1) * double(squeeze(f(patch.S.node)))');
        fun = UserDefinedFunction(fun,d,patch.S.nbnode);
        fun.evaluationAtMultiplePoints = true;
        
        c_patch{k} = FENODEFIELD(H.projection(fun,I));
        
        fun = @(xi) R1_out * (ones(size(xi,1),patch.S.nbnode) + 0.25 * (2* xi(:,3*k) - 1) * double(squeeze(f(patch.S.node)))');
        fun = UserDefinedFunction(fun,d,patch.S.nbnode);
        fun.evaluationAtMultiplePoints = true;
        
        R_patch{k} = FENODEFIELD(H.projection(fun,I));
        
        K_in{k} = K_out;
        c_in{k} = c_out;
        R_in{k} = R1_out;
    end
    
    % Complementary subdomain
    mat_out = MATERIALS();
    mat_out{1} = FOUR_ISOT('k',K_out,'c',c_out,'b',V_out,'r',R1_out);
    mat_out{2} = FOUR_ISOT('k',K_out,'c',c_out,'b',V_out,'r',R2_out);
    mat_out{1} = setnumber(mat_out{1},1);
    mat_out{2} = setnumber(mat_out{2},2);
    glob.S = setmaterial(glob.S,mat_out{1},1);
    glob.S = setmaterial(glob.S,mat_out{2},2:3);
    
    % Patches
    mat_patch = MATERIALS();
    for k=1:n
        mat_patch{k} = FOUR_ISOT('k',K_patch{k},'c',c_patch{k},'b',V_patch{k},'r',R_patch{k});
        mat_patch{k} = setnumber(mat_patch{k},2+k);
        patches.patches{k}.S = setmaterial(patches.patches{k}.S,mat_patch{k});
    end
    
    % Fictitious patches
    mat_in = MATERIALS();
    for k=1:n
        mat_in{k} = FOUR_ISOT('k',K_in{k},'c',c_in{k},'b',V_in{k},'r',R_in{k});
        mat_in{k} = setnumber(mat_in{k},2+k);
        glob.S = setmaterial(glob.S,mat_in{k},getnumgroupelemwithparam(glob.S,'partition',k));
    end
    
    %% Dirichlet boundary conditions
    % Global
    numnode = getnumber(getnode(create_boundary(glob.S)));
    [~,numnode1] = intersect(glob.S,L1);
    [~,numnode2] = intersect(glob.S,L2);
    numnoderest = setdiff(setdiff(numnode,numnode1),numnode2);
    
    glob.S = final(glob.S);
    glob.S = addcl(glob.S,numnode1,'T',1);
    glob.S = addcl(glob.S,numnode2,'T',0);
    permeability = true; % if false, the walls of the canister are considered to be impermeable
    if permeability
        glob.S = addcl(glob.S,numnoderest,'T',0);
    end
    glob.S_out = getfinalmodelpart(glob.S,0);
    % S_in = cell(1,n);
    % for k=1:n
    %     S_in{k} = getfinalmodelpart(glob.S,k);
    % end
    
    % Complementary subdomain
    globOut.S = glob.S_out;
    
    % Patches
    for k=1:n
        patches.patches{k}.S = final(patches.patches{k}.S);
    end
    
    % Interfaces
    interfaces = Interfaces(patches);
    
    %% Initial conditions
    % already taken into account by Dirichlet boundary conditions
    % Global
    % glob.u0 = calc_init_dirichlet(glob.S);
    
    % Complementary subdomain
    % globOut.u0 = calc_init_dirichlet(globOut.S);
    
    % Patches
    % for k=1:n
    %     patches.patches{k}.u0 = calc_init_dirichlet(patches.patches{k}.S);
    % end
    
    %% Time scheme
    t0 = 0;
    t1 = 2;
    nt = 100;
    T = TIMEMODEL(t0,t1,nt);
    
    % N = EULERTIMESOLVER(T,'eulertype','explicit','display',false);
    N = EULERTIMESOLVER(T,'eulertype','implicit','display',false);
    % N = DGTIMESOLVER(T,1,'outputsplit',true,'display',false,'lu',true);
    
    loadFunction = @(N) one(N);
    
    % Global
    glob.timeSolver = N;
    glob.timeOrder = 1;
    
    % Complementary subdomain
    globOut.timeSolver = N;
    globOut.timeOrder = 1;
    
    % Patches
    for k=1:n
        patches.patches{k}.timeSolver = N;
        patches.patches{k}.timeOrder = 1;
    end
    
    %% Mass and stifness matrices and sollicitation vectors
    % Global
    glob.M = calc_mass(glob.S);
    glob.A = calc_rigi(glob.S);
    [~,glob.b0_out] = calc_rigi(glob.S,'selgroup',getnumgroupelemwithparam(glob.S,'partition',0));
    glob.b0_out = -glob.b0_out;
    for k=1:n
        glob.M_in{k} = calc_mass(glob.S,'selgroup',getnumgroupelemwithparam(glob.S,'partition',k));
        glob.A_in{k} = calc_rigi(glob.S,'selgroup',getnumgroupelemwithparam(glob.S,'partition',k));
    end
    glob.b_out = glob.b0_out*loadFunction(glob.timeSolver);
    
    % Complementary subdomain
    globOut.M = calc_mass(globOut.S);
    [globOut.A,globOut.b0] = calc_rigi(globOut.S);
    globOut.b0 = -globOut.b0;
    globOut.b = globOut.b0*loadFunction(globOut.timeSolver);
    
    % Patches
    for k=1:n
        if ~israndom(patches.patches{k}.S)
            patches.patches{k}.M = calc_mass(patches.patches{k}.S);
            [patches.patches{k}.A,patches.patches{k}.b0] = calc_rigi(patches.patches{k}.S);
            patches.patches{k}.b0 = -patches.patches{k}.b0;
            patches.patches{k}.b = patches.patches{k}.b0*loadFunction(patches.patches{k}.timeSolver);
        end
    end
    
    %% Mass matrices
    for k=1:n
        interfaces.interfaces{k}.M = calc_massgeom(interfaces.interfaces{k}.S);
    end
    
    %% Projection operators
    glob.P_out = calcProjection(glob);
    for k=1:n
        interfaces.interfaces{k}.P_glob = calcProjection(glob,interfaces.interfaces{k});
        interfaces.interfaces{k}.P_globOut = interfaces.interfaces{k}.P_glob*glob.P_out';
        interfaces.interfaces{k}.P_patch = calcProjection(patches.patches{k},interfaces.interfaces{k});
    end
    
    %% Parameters for global and local problems
    % Global problem
    glob.increment = false;
    
    % Local problems
    for k=1:n
        patches.patches{k}.changeOfVariable = false;
        patches.patches{k}.increment = false;
    end
    
    %% Stationary solution
    glob_sta = glob;
    glob_sta.timeSolver = [];
    glob_sta.timeOrder = 0;
    glob_sta.b_out = glob.b0_out;
    
    globOut_sta = globOut;
    globOut_sta.timeSolver = [];
    globOut_sta.timeOrder = 0;
    globOut_sta.b = globOut.b0;
    
    patches_sta = patches;
    interfaces_sta = interfaces;
    for k=1:n
        patches_sta.patches{k}.timeSolver = [];
        patches_sta.patches{k}.timeOrder = 0;
        patches_sta.patches{k}.b = patches.patches{k}.b0;
    end
    
    %% Save variables
    save(fullfile(pathname,'problem.mat'),'glob_sta','globOut_sta','patches_sta','interfaces_sta','rv');
    save(fullfile(pathname,'problem_time.mat'),'glob','globOut','patches','interfaces','N','D_patch','Sadv','Sadv_patch','v','v_patch','phi','phi_patch');
else
    load(fullfile(pathname,'problem.mat'),'glob_sta','globOut_sta','patches_sta','interfaces_sta','rv');
    load(fullfile(pathname,'problem_time.mat'),'glob','globOut','patches','interfaces','N','D_patch','Sadv','Sadv_patch','v','v_patch','phi','phi_patch');
end

%% Direct solver
if directSolver
    p = 50;
    bases = cellfun(@(x) PolynomialFunctionalBasis(x,0:p),orthonormalPolynomials(rv),'UniformOutput',false);
    bases = FunctionalBases(bases);
    
    s = AdaptiveSparseTensorAlgorithm();
    s.tolStagnation = 1e-1;
    s.display = true;
    s.displayIterations = true;
    
    ls = LinearModelLearningSquareLoss();
    ls.errorEstimation = true;
    ls.sharedCoefficients = false;
    
    DS = DirectSolver();
    DS.changeOfVariable = false;
    DS.display = true;
    
    % Stationary solution
    s.tol = 1e-6;
    DS.timeSolver = [];
    DS.timeOrder = 0;
    [U_ref,w_ref,lambda_ref,output_ref] = DS.solveRandom(globOut_sta,patches_sta,interfaces_sta,s,bases,ls,rv);
    
    % Transient solution
    s.tol = 1e-3;
    DS.timeSolver = N;
    DS.timeOrder = 1;
    [Ut_ref,wt_ref,lambdat_ref,outputt_ref] = DS.solveRandom(globOut,patches,interfaces,s,bases,ls,rv);
    
    save(fullfile(pathname,'reference_solution.mat'),'U_ref','w_ref','lambda_ref','output_ref','s','bases','ls');
    save(fullfile(pathname,'reference_solution_time.mat'),'Ut_ref','wt_ref','lambdat_ref','outputt_ref','s','bases','ls');
else
    load(fullfile(pathname,'reference_solution.mat'),'U_ref','w_ref','lambda_ref','output_ref','s','bases','ls');
    load(fullfile(pathname,'reference_solution_time.mat'),'Ut_ref','wt_ref','lambdat_ref','outputt_ref','s','bases','ls');
end

%% Outputs
fprintf('\n')
fprintf('Reference stationary solution\n');
fprintf('-----------------------------\n');
dispMultiscaleSolution(U_ref,w_ref,lambda_ref);
fprintf('nb samples = %d\n',output_ref.nbSamples)
fprintf('CV error = %d for U\n',norm(output_ref.CVErrorGlobalSolution))
for k=1:n
    fprintf('         = %d for w{%u}\n',norm(output_ref.CVErrorLocalSolution{k}),k)
    fprintf('         = %d for lambda{%u}\n',norm(output_ref.CVErrorLagrangeMultiplier{k}),k)
end
fprintf('elapsed time = %f s\n',output_ref.time)

fprintf('\n')
fprintf('Reference transient solution\n');
fprintf('----------------------------\n');
dispMultiscaleSolution(Ut_ref,wt_ref,lambdat_ref);
fprintf('nb samples = %d\n',outputt_ref.nbSamples)
fprintf('CV error = %d for U_ref\n',norm(outputt_ref.CVErrorGlobalSolution))
for k=1:n
    fprintf('         = %d for w_ref{%u}\n',norm(outputt_ref.CVErrorLocalSolution{k}),k)
    fprintf('         = %d for lambda_ref{%u}\n',norm(outputt_ref.CVErrorLagrangeMultiplier{k}),k)
end
fprintf('elapsed time = %f s\n',outputt_ref.time)

%% Global-local Iterative solver
if iterativeSolver
    s.tolStagnation = 1e-1;
    s.display = true;
    s.displayIterations = false;
    
    IS = IterativeSolver();
    IS.tolerance = eps;
    IS.relaxation = 'Aitken';
    IS.updateRelaxationParameter = true;
    IS.errorCriterion = 'reference';
    IS.display = true;
    IS.displayIterations = true;
    
    % Stationary solution
    s.tol = 1e-3;
    % s.tol = 10^(-tol-2);
    IS.maxIterations = 20;
    IS.referenceSolution = {U_ref,w_ref,lambda_ref};
    [U,w,lambda,output] = IS.solveRandom(glob_sta,patches_sta,interfaces,s,bases,ls,rv);
    
    % Transient solution
    s.tol = 5e-2;
    % s.tol = 10^(-tol);
    IS.maxIterations = 20;
    IS.referenceSolution = {Ut_ref,wt_ref,lambdat_ref};
    [Ut,wt,lambdat,outputt] = IS.solveRandom(glob,patches,interfaces,s,bases,ls,rv);
    
    save(fullfile(pathname,'solution.mat'),'U','w','lambda','output');
    save(fullfile(pathname,'solution_time.mat'),'Ut','wt','lambdat','outputt');
else
    load(fullfile(pathname,'solution.mat'),'U','w','lambda','output');
    load(fullfile(pathname,'solution_time.mat'),'Ut','wt','lambdat','outputt');
end

%% Outputs
fprintf('\n')
fprintf('Stationary solution\n');
fprintf('-------------------\n');
dispMultiscaleSolution(U,w,lambda);
fprintf('elapsed time = %f s\n',output.totalTime)

fprintf('\n')
fprintf('Transient solution\n');
fprintf('------------------\n');
dispMultiscaleSolution(Ut,wt,lambdat);
fprintf('elapsed time = %f s\n',outputt.totalTime)

%% Display
if displaySolution
    close all
    
    %% Display domains and meshes
    figure('Name','Domain')
    clf
    plot(create_boundary(glob.S));
    hold on
    h1 = plot(glob.S,'selgroup',[1,3+(1:numel(patches))],'FaceColor',getfacecolor(1),'EdgeColor','none');
    h2 = plot(glob.S,'selgroup',2,'FaceColor',[0 0 1],'EdgeColor','none');
    h3 = plot(glob.S,'selgroup',3,'FaceColor',[0 1 0],'EdgeColor','none');
    h4 = plotfacets(glob.S,5,'FaceColor',[0.63 0.13 0.94],'EdgeColor','none');
    h5 = plotfacets(glob.S,17,'FaceColor',[1 0.5 0],'EdgeColor','none');
    hold off
    set(gca,'FontSize',fontsize)
    legend([h1(1),h2(1),h3(1),h4(1),h5(1)],'$\Omega_1$','$\Omega_2$','$\Omega_3$','$\Gamma_D^1$','$\Gamma_D^2$',...
        'Location','NorthEast','Interpreter','latex')
    mysaveas(pathname,'domain',formats,renderer);
    mymatlab2tikz(pathname,'domain.tex');
    
    figure('Name','Advection velocity')
    clf
    plot(phi,Sadv);
    colorbar
    set(gca,'FontSize',fontsize)
    hold on
    ampl = 6;
    quiver(v,Sadv,ampl,'k');
    hold off
    ylim([0,1.7])
    mysaveas(pathname,'advection_velocity_global',formats,renderer);
    
    figure('Name','Advection velocity')
    clf
    Sadv_out = getfinalmodelpart(Sadv,0);
    phi_out = calc_P_free(Sadv,Sadv_out)*phi;
    plot(phi_out,Sadv_out);
    hold on
    for k=1:n
        plot(phi_patch{k},Sadv_patch{k});
    end
    colorbar
    set(gca,'FontSize',fontsize)
    ampl = 6;
    v_out = calc_P(Sadv,Sadv_out)*v;
    quiver(v_out,Sadv_out,ampl,'k');
    for k=1:n
        quiver(v_patch{k},Sadv_patch{k},ampl,'k');
    end
    hold off
    ylim([0,1.7])
    mysaveas(pathname,'advection_velocity_global_patches',formats,renderer);
    
    % plotDomain(glob.S,D_patch);
    % mysaveas(pathname,'domain_global_patches',formats,renderer);
    % mymatlab2tikz(pathname,'domain_global_patches.tex');
    
    numbers = getnumber(patches);
    figure('Name',['Complement subdomain and patches #' num2str([numbers{:}])])
    clf
    plot(create_boundary(glob.S));
    hold on
    h1 = plot(glob.S,'selgroup',1:3,'FaceColor',getfacecolor(1),'EdgeColor','none');
    h2 = plot(glob.S,'selgroup',2,'FaceColor',[0 0 1],'EdgeColor','none');
    h3 = plot(glob.S,'selgroup',3,'FaceColor',[0 1 0],'EdgeColor','none');
    h4 = plotfacets(glob.S,5,'FaceColor',[0.63 0.13 0.94],'EdgeColor','none');
    h5 = plotfacets(glob.S,17,'FaceColor',[1 0.5 0],'EdgeColor','none');
    h_patch = cell(1,n);
    leg_patch = cell(1,n);
    for k=1:n
        plot(create_boundary(patches.patches{k}.S));
        h_patch{k} = plot(patches.patches{k}.S,'FaceColor',getfacecolor(1+k),'EdgeColor','none');
        leg_patch{k} = ['$\Lambda_{' num2str(k) '}$'];
    end
    hold off
    set(gca,'FontSize',fontsize)
    legend([h1(1),h2(1),h3(1),h4(1),h5(1),h_patch{:}],...
        '$\Omega_1 \setminus \Lambda$','$\Omega_2$','$\Omega_3$','$\Gamma_D^1$','$\Gamma_D^2$',leg_patch{:},...
        'Location','NorthEast','Interpreter','latex')
    mysaveas(pathname,'domain_global_patches',formats,renderer);
    mymatlab2tikz(pathname,'domain_global_patches.tex');
    
    % plotPartition(glob,'legend',false);
    % mysaveas(pathname,'mesh_partition',formats,renderer);
    
    % plotModel(glob,patches,'legend',false);
    % mysaveas(pathname,'mesh_global_patches',formats,renderer);
    
    figure('Name',['Meshes of complementary subdomain and patches #' num2str([numbers{:}])])
    clf
    h1 = plot(glob.S,'selgroup',1:3,'FaceColor',getfacecolor(1));
    hold on
    h2 = plot(glob.S,'selgroup',2,'FaceColor',[0 0 1]);
    h3 = plot(glob.S,'selgroup',3,'FaceColor',[0 1 0]);
    h4 = plotfacets(glob.S,5,'FaceColor',[0.63 0.13 0.94]);
    h5 = plotfacets(glob.S,17,'FaceColor',[1 0.5 0]);
    h_patch = cell(1,n);
    leg_patch = cell(1,n);
    for k=1:n
        h_patch{k} = plot(patches.patches{k}.S,'FaceColor',getfacecolor(1+k));
        leg_patch{k} = ['$\Lambda_{' num2str(k) '}$'];
    end
    hold off
    set(gca,'FontSize',fontsize)
    legend([h1(1),h2(1),h3(1),h4(1),h5(1),h_patch{:}],...
        '$\Omega_1 \setminus \Lambda$','$\Omega_2$','$\Omega_3$','$\Gamma_D^1$','$\Gamma_D^2$',leg_patch{:},...
        'Location','NorthEast','Interpreter','latex')
    mysaveas(pathname,'mesh_global_patches',formats,renderer);
    
    %% Display evolutions of error indicator, stagnation indicator, CPU time w.r.t. number of iterations for stationary solution
    plotError(output);
    mysaveas(pathname,'error_stationary','fig');
    mymatlab2tikz(pathname,'error_stationary.tex');
    
    plotStagnation(output);
    mysaveas(pathname,'stagnation_stationary','fig');
    mymatlab2tikz(pathname,'stagnation_stationary.tex');
    
    plotErrorGlobalSolution(output);
    mysaveas(pathname,'error_global_solution_stationary','fig');
    mymatlab2tikz(pathname,'error_global_solution_stationary.tex');
    
    plotStagnationGlobalSolution(output);
    mysaveas(pathname,'stagnation_global_solution_stationary','fig');
    mymatlab2tikz(pathname,'stagnation_global_solution_stationary.tex');
    
    plotCPUTime(output,'legend',false);
    mysaveas(pathname,'cpu_time_stationary','fig');
    mymatlab2tikz(pathname,'cpu_time_stationary.tex');
    
    plotRelaxationParameter(output,'legend',false);
    mysaveas(pathname,'relaxation_parameter_stationary','fig');
    mymatlab2tikz(pathname,'relaxation_parameter_stationary.tex');
    
    plotNbSamples(output);
    mysaveas(pathname,'nb_samples_stationary','fig');
    mymatlab2tikz(pathname,'nb_samples_stationary.tex');
    
    plotDimStochasticBasis(output);
    mysaveas(pathname,'dim_stochastic_basis_stationary','fig');
    mymatlab2tikz(pathname,'dim_stochastic_basis_stationary.tex');
    
    plotCVError(output);
    mysaveas(pathname,'cv_error_stationary','fig');
    mymatlab2tikz(pathname,'cv_error_stationary.tex');
    
    %% Display statistical outputs for stationary solutions
    % plotStatsAllSolutions(glob,patches,interfaces,U,w,lambda);
    % plotStatsAllSolutions(glob,patches,interfaces,U,w,lambda,'surface',true);
    
    plotMeanGlobalSolution(glob,U);
    mysaveas(pathname,'mean_global_solution',formats,renderer);
    plotMeanGlobalSolution(glob,U,'surface',true);
    mysaveas(pathname,'mean_global_solution_surface',formats,renderer);
    
    % plotMeanLocalSolution(patches,w);
    % mysaveas(pathname,'mean_local_solution',formats,renderer);
    % plotMeanLocalSolution(patches,w,'surface',true);
    % mysaveas(pathname,'mean_local_solution_surface',formats,renderer);
    
    % plotMeanLagrangeMultiplier(interfaces,lambda);
    % mysaveas(pathname,'mean_Lagrange_multiplier',formats,renderer);
    % plotMeanLagrangeMultiplier(interfaces,lambda,'surface',true);
    % mysaveas(pathname,'mean_Lagrange_multiplier_surface',formats,renderer);
    
    plotMeanMultiscaleSolution(glob,patches,interfaces,U,w);
    mysaveas(pathname,'mean_multiscale_solution',formats,renderer);
    plotMeanMultiscaleSolution(glob,patches,interfaces,U,w,'surface',true);
    mysaveas(pathname,'mean_multiscale_solution_surface',formats,renderer);
    
    plotMeanGlobalLocalSolution(glob,patches,interfaces,U,w);
    mysaveas(pathname,'mean_global_local_solution',formats,renderer);
    plotMeanGlobalLocalSolution(glob,patches,interfaces,U,w,'surface',true);
    mysaveas(pathname,'mean_global_local_solution_surface',formats,renderer);
    
    plotVarianceGlobalSolution(glob,U);
    mysaveas(pathname,'var_global_solution',formats,renderer);
    plotVarianceGlobalSolution(glob,U,'surface',true);
    mysaveas(pathname,'var_global_solution_surface',formats,renderer);
    
    % plotVarianceLocalSolution(patches,w);
    % mysaveas(pathname,'var_local_solution',formats,renderer);
    % plotVarianceLocalSolution(patches,w,'surface',true);
    % mysaveas(pathname,'var_local_solution_surface',formats,renderer);
    
    % plotVarianceLagrangeMultiplier(interfaces,lambda);
    % mysaveas(pathname,'var_Lagrange_multiplier',formats,renderer);
    % plotVarianceLagrangeMultiplier(interfaces,lambda,'surface',true);
    % mysaveas(pathname,'var_Lagrange_multiplier_surface',formats,renderer);
    
    plotVarianceMultiscaleSolution(glob,patches,interfaces,U,w);
    mysaveas(pathname,'var_multiscale_solution',formats,renderer);
    plotVarianceMultiscaleSolution(glob,patches,interfaces,U,w,'surface',true);
    mysaveas(pathname,'var_multiscale_solution_surface',formats,renderer);
    
    plotVarianceGlobalLocalSolution(glob,patches,interfaces,U,w);
    mysaveas(pathname,'var_global_local_solution',formats,renderer);
    plotVarianceGlobalLocalSolution(glob,patches,interfaces,U,w,'surface',true);
    mysaveas(pathname,'var_global_local_solution_surface',formats,renderer);
    
    plotStdGlobalSolution(glob,U);
    mysaveas(pathname,'std_global_solution',formats,renderer);
    plotStdGlobalSolution(glob,U,'surface',true);
    mysaveas(pathname,'std_global_solution_surface',formats,renderer);
    
    % plotStdLocalSolution(patches,w);
    % mysaveas(pathname,'std_local_solution',formats,renderer);
    % plotStdLocalSolution(patches,w,'surface',true);
    % mysaveas(pathname,'std_local_solution_surface',formats,renderer);
    
    % plotStdLagrangeMultiplier(interfaces,lambda);
    % mysaveas(pathname,'std_Lagrange_multiplier',formats,renderer);
    % plotStdLagrangeMultiplier(interfaces,lambda,'surface',true);
    % mysaveas(pathname,'std_Lagrange_multiplier_surface',formats,renderer);
    
    plotStdMultiscaleSolution(glob,patches,interfaces,U,w);
    mysaveas(pathname,'std_multiscale_solution',formats,renderer);
    plotStdMultiscaleSolution(glob,patches,interfaces,U,w,'surface',true);
    mysaveas(pathname,'std_multiscale_solution_surface',formats,renderer);
    
    plotStdGlobalLocalSolution(glob,patches,interfaces,U,w);
    mysaveas(pathname,'std_global_local_solution',formats,renderer);
    plotStdGlobalLocalSolution(glob,patches,interfaces,U,w,'surface',true);
    mysaveas(pathname,'std_global_local_solution_surface',formats,renderer);
    
    d = ndims(U.basis);
    for i=1:d
        % plotSobolIndicesMultiscaleSolution(glob,patches,interfaces,U,w,i);
        % mysaveas(pathname,['sobol_indices_multiscale_solution_var_' num2str(i)],formats,renderer);
        % plotSobolIndicesMultiscaleSolution(glob,patches,interfaces,U,w,i,'surface',true);
        % mysaveas(pathname,['sobol_indices_multiscale_solution_var_' num2str(i) '_surface'],formats,renderer);
        
        plotSensitivityIndicesMultiscaleSolution(glob,patches,interfaces,U,w,i);
        mysaveas(pathname,['sensitivity_indices_multiscale_solution_var_' num2str(i)],formats,renderer);
        plotSensitivityIndicesMultiscaleSolution(glob,patches,interfaces,U,w,i,'surface',true);
        mysaveas(pathname,['sensitivity_indices_multiscale_solution_var_' num2str(i) '_surface'],formats,renderer);
    end
    
    close all
    
    %% Display evolutions of error indicator, stagnation indicator, CPU time w.r.t. number of iterations for transient solutions
    plotError(outputt);
    mysaveas(pathname,'error_transient','fig');
    mymatlab2tikz(pathname,'error_transient.tex');
    
    plotStagnation(outputt);
    mysaveas(pathname,'stagnation_transient','fig');
    mymatlab2tikz(pathname,'stagnation_transient.tex');
    
    plotErrorGlobalSolution(outputt);
    mysaveas(pathname,'error_global_solution_transient','fig');
    mymatlab2tikz(pathname,'error_global_solution_transient.tex');
    
    plotStagnationGlobalSolution(outputt);
    mysaveas(pathname,'stagnation_global_solution_transient','fig');
    mymatlab2tikz(pathname,'stagnation_global_solution_transient.tex');
    
    plotCPUTime(outputt,'legend',false);
    mysaveas(pathname,'cpu_time_transient','fig');
    mymatlab2tikz(pathname,'cpu_time_transient.tex');
    
    plotRelaxationParameter(outputt,'legend',false);
    mysaveas(pathname,'relaxation_parameter_transient','fig');
    mymatlab2tikz(pathname,'relaxation_parameter_transient.tex');
    
    plotNbSamples(outputt);
    mysaveas(pathname,'nb_samples_transient','fig');
    mymatlab2tikz(pathname,'nb_samples_transient.tex');
    
    plotDimStochasticBasis(outputt);
    mysaveas(pathname,'dim_stochastic_basis_transient','fig');
    mymatlab2tikz(pathname,'dim_stochastic_basis_transient.tex');
    
    plotCVError(outputt);
    mysaveas(pathname,'cv_error_transient','fig');
    mymatlab2tikz(pathname,'cv_error_transient.tex');
    
    %% Display evolutions of statistical outputs for transient solutions
    T = gettimemodel(glob.timeSolver);
    
    evolMeanGlobalSolution(glob,T,Ut,'filename','mean_global_solution','pathname',pathname);
    evolMeanGlobalSolution(glob,T,Ut,'surface',true,'filename','mean_global_solution_surface','pathname',pathname);
    
    % evolMeanLocalSolution(patches,T,wt,'filename','mean_local_solution','pathname',pathname);
    % evolMeanLocalSolution(patches,T,wt,'surface',true,'filename','mean_local_solution_surface','pathname',pathname);
    
    % evolMeanLagrangeMultiplier(interfaces,T,lambdat,'filename','mean_Lagrange_multiplier','pathname',pathname);
    % evolMeanLagrangeMultiplier(interfaces,T,lambdat,'surface',true,'filename','mean_Lagrange_multiplier_surface','pathname',pathname);
    
    evolMeanMultiscaleSolution(glob,patches,interfaces,T,Ut,wt,'filename','mean_multiscale_solution','pathname',pathname);
    evolMeanMultiscaleSolution(glob,patches,interfaces,T,Ut,wt,'surface',true,'filename','mean_multiscale_solution_surface','pathname',pathname);
    
    evolMeanGlobalLocalSolution(glob,patches,interfaces,T,Ut,wt,'filename','mean_global_local_solution','pathname',pathname);
    evolMeanGlobalLocalSolution(glob,patches,interfaces,T,Ut,wt,'surface',true,'filename','mean_global_local_solution_surface','pathname',pathname);
    
    evolVarianceGlobalSolution(glob,T,Ut,'filename','var_global_solution','pathname',pathname);
    evolVarianceGlobalSolution(glob,T,Ut,'surface',true,'filename','var_global_solution_surface','pathname',pathname);
    
    % evolVarianceLocalSolution(patches,T,wt,'filename','var_local_solution','pathname',pathname);
    % evolVarianceLocalSolution(patches,T,wt,'surface',true,'filename','var_local_solution_surface','pathname',pathname);
    
    % evolVarianceLagrangeMultiplier(interfaces,T,lambdat,'filename','var_Lagrange_multiplier','pathname',pathname);
    % evolVarianceLagrangeMultiplier(interfaces,T,lambdat,'surface',true,'filename','var_Lagrange_multiplier_surface','pathname',pathname);
    
    evolVarianceMultiscaleSolution(glob,patches,interfaces,T,Ut,wt,'filename','var_multiscale_solution','pathname',pathname);
    evolVarianceMultiscaleSolution(glob,patches,interfaces,T,Ut,wt,'surface',true,'filename','var_multiscale_solution_surface','pathname',pathname);
    
    evolVarianceGlobalLocalSolution(glob,patches,interfaces,T,Ut,wt,'filename','var_global_local_solution','pathname',pathname);
    evolVarianceGlobalLocalSolution(glob,patches,interfaces,T,Ut,wt,'surface',true,'filename','var_global_local_solution_surface','pathname',pathname);
    
    evolStdGlobalSolution(glob,T,Ut,'filename','std_global_solution','pathname',pathname);
    evolStdGlobalSolution(glob,T,Ut,'surface',true,'filename','std_global_solution_surface','pathname',pathname);
    
    % evolStdLocalSolution(patches,T,wt,'filename','std_local_solution','pathname',pathname);
    % evolStdLocalSolution(patches,T,wt,'surface',true,'filename','std_local_solution_surface','pathname',pathname);
    
    % evolStdLagrangeMultiplier(interfaces,T,lambdat,'filename','std_Lagrange_multiplier','pathname',pathname);
    % evolStdLagrangeMultiplier(interfaces,T,lambdat,'surface',true,'filename','std_Lagrange_multiplier_surface','pathname',pathname);
    
    evolStdMultiscaleSolution(glob,patches,interfaces,T,Ut,wt,'filename','std_multiscale_solution','pathname',pathname);
    evolStdMultiscaleSolution(glob,patches,interfaces,T,Ut,wt,'surface',true,'filename','std_multiscale_solution_surface','pathname',pathname);
    
    evolStdGlobalLocalSolution(glob,patches,interfaces,T,Ut,wt,'filename','std_global_local_solution','pathname',pathname);
    evolStdGlobalLocalSolution(glob,patches,interfaces,T,Ut,wt,'surface',true,'filename','std_global_local_solution_surface','pathname',pathname);
    
    d = ndims(Ut.basis);
    for i=1:d
        % evolSobolIndicesMultiscaleSolution(glob,patches,interfaces,T,Ut,wt,i,'filename',['sobol_indices_multiscale_solution_var_' num2str(i)],'pathname',pathname);
        % evolSobolIndicesMultiscaleSolution(glob,patches,interfaces,T,Ut,wt,i,'surface',true,'filename',['sobol_indices_multiscale_solution_var_' num2str(i) '_surface'],'pathname',pathname);
        
        evolSensitivityIndicesMultiscaleSolution(glob,patches,interfaces,T,Ut,wt,i,'filename',['sensitivity_indices_multiscale_solution_var_' num2str(i)],'pathname',pathname);
        evolSensitivityIndicesMultiscaleSolution(glob,patches,interfaces,T,Ut,wt,i,'surface',true,'filename',['sensitivity_indices_multiscale_solution_var_' num2str(i) '_surface'],'pathname',pathname);
    end
    
    %% Display statistical outputs for transient solutions at different instants
%     [t,rep] = gettevol(N);
%     for k=1:floor(length(rep)/4):length(rep)
%         close all
%         Uk_data = Ut.data(:,:,rep(k));
%         wk_data = cellfun(@(x) x.data(:,:,rep(k)),wt,'UniformOutput',false);
%         lambdak_data = cellfun(@(x) x.data(:,:,rep(k)),lambdat,'UniformOutput',false);
%         Uk = FunctionalBasisArray(Uk_data,Ut.basis,Ut.sz(1));
%         wk = cellfun(@(data,x) FunctionalBasisArray(data,x.basis,x.sz(1)),wk_data,wt,'UniformOutput',false);
%         lambdak = cellfun(@(data,x) FunctionalBasisArray(data,x.basis,x.sz(1)),lambdak_data,lambdat,'UniformOutput',false);
%         
%         % plotStatsAllSolutions(glob,patches,interfaces,Uk,wk,lambdak);
%         % plotStatsAllSolutions(glob,patches,interfaces,Uk,wk,lambdak,'surface',true);
%         
%         plotMeanGlobalSolution(glob,Uk);
%         mysaveas(pathname,['mean_global_solution_t' num2str(k-1)],formats,renderer);
%         plotMeanGlobalSolution(glob,Uk,'surface',true);
%         mysaveas(pathname,['mean_global_solution_t' num2str(k-1) '_surface'],formats,renderer);
%         
%         % plotMeanLocalSolution(patches,wk);
%         % mysaveas(pathname,['mean_local_solution_t' num2str(k-1)],formats,renderer);
%         % plotMeanLocalSolution(patches,wk,'surface',true);
%         % mysaveas(pathname,['mean_local_solution_t' num2str(k-1) '_surface'],formats,renderer);
%         
%         % plotMeanLagrangeMultiplier(interfaces,lambdak);
%         % mysaveas(pathname,['mean_Lagrange_multiplier_t' num2str(k-1)],formats,renderer);
%         % plotMeanLagrangeMultiplier(interfaces,lambdak,'surface',true);
%         % mysaveas(pathname,['mean_Lagrange_multiplier_t' num2str(k-1) '_surface'],formats,renderer);
%         
%         plotMeanMultiscaleSolution(glob,patches,interfaces,Uk,wk);
%         mysaveas(pathname,['mean_multiscale_solution_t' num2str(k-1)],formats,renderer);
%         plotMeanMultiscaleSolution(glob,patches,interfaces,Uk,wk,'surface',true);
%         mysaveas(pathname,['mean_multiscale_solution_t' num2str(k-1) '_surface'],formats,renderer);
%         
%         plotMeanGlobalLocalSolution(glob,patches,interfaces,Uk,wk);
%         mysaveas(pathname,['mean_global_local_solution_t' num2str(k-1)],formats,renderer);
%         plotMeanGlobalLocalSolution(glob,patches,interfaces,Uk,wk,'surface',true);
%         mysaveas(pathname,['mean_global_local_solution_t' num2str(k-1) '_surface'],formats,renderer);
%         
%         plotVarianceGlobalSolution(glob,Uk);
%         mysaveas(pathname,['var_global_solution_t' num2str(k-1)],formats,renderer);
%         plotVarianceGlobalSolution(glob,Uk,'surface',true);
%         mysaveas(pathname,['var_global_solution_t' num2str(k-1) '_surface'],formats,renderer);
%         
%         % plotVarianceLocalSolution(patches,wk);
%         % mysaveas(pathname,['var_local_solution_t' num2str(k-1)],formats,renderer);
%         % plotVarianceLocalSolution(patches,wk,'surface',true);
%         % mysaveas(pathname,['var_local_solution_t' num2str(k-1) '_surface'],formats,renderer);
%         
%         % plotVarianceLagrangeMultiplier(interfaces,lambdak);
%         % mysaveas(pathname,['var_Lagrange_multiplier_t' num2str(k-1)],formats,renderer);
%         % plotVarianceLagrangeMultiplier(interfaces,lambdak,'surface',true);
%         % mysaveas(pathname,['var_Lagrange_multiplier_t' num2str(k-1) '_surface'],formats,renderer);
%         
%         plotVarianceMultiscaleSolution(glob,patches,interfaces,Uk,wk);
%         mysaveas(pathname,['var_multiscale_solution_t' num2str(k-1)],formats,renderer);
%         plotVarianceMultiscaleSolution(glob,patches,interfaces,Uk,wk,'surface',true);
%         mysaveas(pathname,['var_multiscale_solution_t' num2str(k-1) '_surface'],formats,renderer);
%         
%         plotVarianceGlobalLocalSolution(glob,patches,interfaces,Uk,wk);
%         mysaveas(pathname,['var_global_local_solution_t' num2str(k-1)],formats,renderer);
%         plotVarianceGlobalLocalSolution(glob,patches,interfaces,Uk,wk,'surface',true);
%         mysaveas(pathname,['var_global_local_solution_t' num2str(k-1) '_surface'],formats,renderer);
%         
%         plotStdGlobalSolution(glob,Uk);
%         mysaveas(pathname,['std_global_solution_t' num2str(k-1)],formats,renderer);
%         plotStdGlobalSolution(glob,Uk,'surface',true);
%         mysaveas(pathname,['std_global_solution_t' num2str(k-1) '_surface'],formats,renderer);
%         
%         % plotStdLocalSolution(patches,wk);
%         % mysaveas(pathname,['std_local_solution_t' num2str(k-1)],formats,renderer);
%         % plotStdLocalSolution(patches,wk,'surface',true);
%         % mysaveas(pathname,['std_local_solution_t' num2str(k-1) '_surface'],formats,renderer);
%         
%         % plotStdLagrangeMultiplier(interfaces,lambdak);
%         % mysaveas(pathname,['std_Lagrange_multiplier_t' num2str(k-1)],formats,renderer);
%         % plotStdLagrangeMultiplier(interfaces,lambdak,'surface',true);
%         % mysaveas(pathname,['std_Lagrange_multiplier_t' num2str(k-1) '_surface'],formats,renderer);
%         
%         plotStdMultiscaleSolution(glob,patches,interfaces,Uk,wk);
%         mysaveas(pathname,['std_multiscale_solution_t' num2str(k-1)],formats,renderer);
%         plotStdMultiscaleSolution(glob,patches,interfaces,Uk,wk,'surface',true);
%         mysaveas(pathname,['std_multiscale_solution_t' num2str(k-1) '_surface'],formats,renderer);
%         
%         plotStdGlobalLocalSolution(glob,patches,interfaces,Uk,wk);
%         mysaveas(pathname,['std_global_local_solution_t' num2str(k-1)],formats,renderer);
%         plotStdGlobalLocalSolution(glob,patches,interfaces,Uk,wk,'surface',true);
%         mysaveas(pathname,['std_global_local_solution_t' num2str(k-1) '_surface'],formats,renderer);
%         
%         d = ndims(Uk.basis);
%         for i=1:d
%             % plotSobolIndicesMultiscaleSolution(glob,patches,interfaces,Uk,wk,i);
%             % mysaveas(pathname,['sobol_indices_multiscale_solution_var_' num2str(i) '_t' num2str(k-1)],formats,renderer);
%             % plotSobolIndicesMultiscaleSolution(glob,patches,interfaces,Uk,wk,i,'surface',true);
%             % mysaveas(pathname,['sobol_indices_multiscale_solution_var_' num2str(i) '_t' num2str(k-1) '_surface'],formats,renderer);
%             
%             plotSensitivityIndicesMultiscaleSolution(glob,patches,interfaces,Uk,wk,i);
%             mysaveas(pathname,['sensitivity_indices_multiscale_solution_var_' num2str(i) '_t' num2str(k-1)],formats,renderer);
%             plotSensitivityIndicesMultiscaleSolution(glob,patches,interfaces,Uk,wk,i,'surface',true);
%             mysaveas(pathname,['sensitivity_indices_multiscale_solution_var_' num2str(i) '_t' num2str(k-1) '_surface'],formats,renderer);
%         end
%     end
    
    %% Display quantity of interest
    % boutput: mean of concentration of pollutant captured by the trap domain
    %          (group #2 in mesh) as a function of time
    % Ioutput: mean of total concentration of pollutant captured by the trap domain
    %          (group #2 in mesh) along the complete time evolution,
    %          corresponding to all the pollutant that the actual filter
    %          (group #1 in mesh) is not able to retain
    foutput = bodyload(keepgroupelem(glob.S,2),[],'QN',1,'nofree');
    foutput_ref = bodyload(keepgroupelem(globOut.S,2),[],'QN',1,'nofree');
    
    sz = [getnbddlfree(glob.S),getnbtimedof(T)];
    sz_ref = [getnbddlfree(globOut.S),getnbtimedof(T)];
    
    mean_Ut = mean(Ut);
    mean_Ut = reshape(mean_Ut,sz);
    mean_Ut = TIMEMATRIX(mean_Ut,T);
    mean_Ut = unfreevector(glob.S,mean_Ut);
    mean_Ut_ref = mean(Ut_ref);
    mean_Ut_ref = reshape(mean_Ut_ref,sz_ref);
    mean_Ut_ref = TIMEMATRIX(mean_Ut_ref,T);
    mean_Ut_ref = unfreevector(globOut.S,mean_Ut_ref);
    
    std_Ut = std(Ut);
    std_Ut = reshape(std_Ut,sz);
    std_Ut = TIMEMATRIX(std_Ut,T);
    std_Ut = unfreevector(glob.S,std_Ut)-calc_init_dirichlet(glob.S)*one(T);
    std_Ut_ref = std(Ut_ref);
    std_Ut_ref = reshape(std_Ut_ref,sz_ref);
    std_Ut_ref = TIMEMATRIX(std_Ut_ref,T);
    std_Ut_ref = unfreevector(globOut.S,std_Ut_ref)-calc_init_dirichlet(globOut.S)*one(T);
    
    mean_boutput = foutput'*mean_Ut;
    mean_boutput_ref = foutput_ref'*mean_Ut_ref;
    std_boutput = foutput'*std_Ut;
    std_boutput_ref = foutput_ref'*std_Ut_ref;
    
    probs = [0.025 0.975];
    nbSamples = 100;
    x = random(rv,nbSamples);
    samples_Ut = Ut(x);
    samples_Ut_ref = Ut_ref(x);
    samples_boutput = zeros(nbSamples,Ut.sz(2));
    samples_boutput_ref = zeros(nbSamples,Ut_ref.sz(2));
    for i=1:nbSamples
        sample_Ut = reshape(samples_Ut(i,:,:),Ut.sz);
        sample_Ut = TIMEMATRIX(sample_Ut,T);
        sample_Ut = unfreevector(glob.S,sample_Ut)-calc_init_dirichlet(glob.S)*one(T);
        sample_Ut_ref = reshape(samples_Ut_ref(i,:,:),Ut_ref.sz);
        sample_Ut_ref = TIMEMATRIX(sample_Ut_ref,T);
        sample_Ut_ref = unfreevector(globOut.S,sample_Ut_ref)-calc_init_dirichlet(globOut.S)*one(T);
        samples_boutput(i,:) = foutput'*sample_Ut;
        samples_boutput_ref(i,:) = foutput_ref'*sample_Ut_ref;
    end
    ci_boutput = quantile(samples_boutput,probs);
    ci_boutput_ref = quantile(samples_boutput_ref,probs);
    
    figure('Name','Quantity of interest')
    clf
    [t,rep] = gettplot(T);
    ciplot(ci_boutput(1,:),ci_boutput(2,:),t,'b');
    hold on
    ciplot(ci_boutput_ref(1,:),ci_boutput_ref(2,:),t,'r');
    alpha(0.2)
    plot(mean_boutput,'-b','LineWidth',1);
    plot(mean_boutput_ref,'-r','LineWidth',1);
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('Time [s]')
    ylabel('Concentration of pollutant')
    legend({['$' num2str((probs(2)-probs(1))*100) '\%$ confidence interval - multiscale'],...
        ['$' num2str((probs(2)-probs(1))*100) '\%$ confidence interval - monoscale'],...
        'mean value - multiscale','mean value - monoscale'},'Interpreter','latex')
    mysaveas(pathname,'quantity_of_interest',formats,renderer);
    mymatlab2tikz(pathname,'quantity_of_interest.tex');
    
    mean_Ioutput = integrate(mean_boutput);
    mean_Ioutput_ref = integrate(mean_boutput_ref);
    mean_errOutput = norm(mean_Ioutput-mean_Ioutput_ref)/mean_Ioutput_ref;
    std_Ioutput = integrate(std_boutput);
    std_Ioutput_ref = integrate(std_boutput_ref);
    std_errOutput = norm(std_Ioutput-std_Ioutput_ref)/std_Ioutput_ref;
    lowerci_Ioutput = integrate(TIMEMATRIX(ci_boutput(1,:),T));
    lowerci_Ioutput_ref = integrate(TIMEMATRIX(ci_boutput_ref(1,:),T));
    lowerci_errOutput = norm(lowerci_Ioutput-lowerci_Ioutput_ref)/lowerci_Ioutput_ref;
    upperci_Ioutput = integrate(TIMEMATRIX(ci_boutput(2,:),T));
    upperci_Ioutput_ref = integrate(TIMEMATRIX(ci_boutput_ref(2,:),T));
    upperci_errOutput = norm(upperci_Ioutput-upperci_Ioutput_ref)/upperci_Ioutput_ref;
    fprintf('mean of quantity of interest           = %e\n',mean_Ioutput);
    fprintf('mean of reference quantity of interest = %e\n',mean_Ioutput_ref);
    fprintf('error in mean of quantity of interest  = %e\n',mean_errOutput);
    fprintf('std of quantity of interest            = %e\n',std_Ioutput);
    fprintf('std of reference quantity of interest  = %e\n',std_Ioutput_ref);
    fprintf('error in std of quantity of interest   = %e\n',std_errOutput);
    fprintf('%d%% confidence interval of quantity of interest            = [%e,%e]\n',(probs(2)-probs(1))*100,lowerci_Ioutput,upperci_Ioutput);
    fprintf('%d%% confidence interval of reference quantity of interest  = [%e,%e]\n',(probs(2)-probs(1))*100,lowerci_Ioutput_ref,upperci_Ioutput_ref);
    fprintf('error in %d%% confidence interval of quantity of interest   = [%e,%e]\n',(probs(2)-probs(1))*100,lowerci_errOutput,upperci_errOutput);
end

% end

myparallel('stop');
