%% Multiscale deterministic transient linear advection-diffusion-reaction problem %%
%%--------------------------------------------------------------------------------%%

% clc
clearvars
close all
% myparallel('start');

%% Input data
setProblem = true;
directSolver = true;
iterativeSolver = true;
displaySolution = true;

n = 3; % number of patches
filename = ['dynLinElas' num2str(n) 'Patches'];
pathname = fullfile(getfemobjectoptions('path'),'MULTISCALE',...
    'examples','multiscaleDet',filename);
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
    
    D = DOMAIN(2,[0.0,0.0],[2.0,0.4]);
    
    elemtype = 'TRI3';
    % elemtype = 'QUA4';
    % option = 'DEFO'; % plane strain
    option = 'CONT'; % plane stress
    nbelem = [60,12];
    glob.S = build_model(D,'nbelem',nbelem,'elemtype',elemtype,'option',option);
%     cl = 0.03;
%     glob.S = build_model(D,'cl',cl,'elemtype',elemtype,'option',option,'filename',fullfile(pathname,'gmsh_domain'));
    
    % Patches
    patches = Patches(n);
    
    D_patch = cell(1,n);
    D_patch{1} = DOMAIN(2,[1.7,0.1],[1.9,0.3]);
    D_patch{2} = DOMAIN(2,[0.9,0.1],[1.1,0.3]);
    D_patch{3} = DOMAIN(2,[0.1,0.1],[0.3,0.3]);
    
    nbelem_patch = [18,18];
    for k=1:n
        patches.patches{k}.S = build_model(D_patch{k},'nbelem',nbelem_patch);
    end
    % cl_patch = 0.005;
    % for k=1:n
    %     patches.patches{k}.S = build_model(D_patch{k},'cl',cl_patch,'filename',fullfile(pathname,['gmsh_patch_' num2str(k)]));
    % end
    
    % Partition of global mesh
    glob = partition(glob,D_patch);
    
    %% Materials
    % Poisson ratio
    NU = 0.3;
    % Thickness
    DIM3 = 1;
    % Density
    RHO_out = 1;
    RHO_patch = cell(1,n);
    RHO_in = cell(1,n);
    % Young modulus
    E_out = 1;
    E_patch = cell(1,n);
    E_in = cell(1,n);
    for k=1:n
        RHO_patch{k} = RHO_out;
        E_patch{k} = E_out;
        RHO_in{k} = RHO_out;
        E_in{k} = E_out;
    end
    
    % Complementary subdomain
    mat_out = ELAS_ISOT('E',E_out,'NU',NU,'RHO',RHO_out,'DIM3',DIM3);
    mat_out = setnumber(mat_out,0);
    glob.S = setmaterial(glob.S,mat_out,getnumgroupelemwithparam(glob.S,'partition',0));
    
    % Patches
    mat_patch = MATERIALS();
    for k=1:n
        mat_patch{k} = ELAS_ISOT('E',E_patch{k},'NU',NU,'RHO',RHO_patch{k},'DIM3',DIM3);
        mat_patch{k} = setnumber(mat_patch{k},k);
        patches.patches{k}.S = setmaterial(patches.patches{k}.S,mat_patch{k});
    end
    
    % Fictitious patches
    mat_in = MATERIALS();
    for k=1:n
        mat_in{k} = ELAS_ISOT('E',E_in{k},'NU',NU,'RHO',RHO_in{k},'DIM3',DIM3);
        mat_in{k} = setnumber(mat_in{k},k);
        glob.S = setmaterial(glob.S,mat_in{k},getnumgroupelemwithparam(glob.S,'partition',k));
    end
    
    %% Dirichlet boundary conditions
    % Global
    L1 = LIGNE(getvertex(D,1),getvertex(D,4));
    L2 = LIGNE(getvertex(D,2),getvertex(D,3));
    
    glob.S = final(glob.S);
    glob.S = addcl(glob.S,L1);
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
    % Global
    glob.u0 = zeros(getnbddlfree(glob.S),1);
    glob.v0 = zeros(getnbddlfree(glob.S),1);
    
    % Complementary subdomain
    globOut.u0 = zeros(getnbddlfree(globOut.S),1);
    globOut.v0 = zeros(getnbddlfree(globOut.S),1);
    
    % Patches
    for k=1:n
        patches.patches{k}.u0 = zeros(getnbddlfree(patches.patches{k}.S),1);
        patches.patches{k}.v0 = zeros(getnbddlfree(patches.patches{k}.S),1);
    end
    
    %% Time scheme
    t0 = 0;
    t1 = 2;
    nt = 50;
    T = TIMEMODEL(t0,t1,nt);
    
    N = NEWMARKSOLVER(T,'alpha',0.05,'display',false);
    % N = DGTIMESOLVER(T,1,'outputsplit',true,'display',false,'lu',true);
    
    tc = get(T,'t1')/6;
    loadFunction = @(N) rampe(N,t0,tc);
    % loadFunction = @(N) dirac(N,t0,tc);
    % loadFunction = @(N) one(N);
    
    % Global
    glob.timeSolver = N;
    glob.timeOrder = 2;
    
    % Complementary subdomain
    globOut.timeSolver = N;
    globOut.timeOrder = 2;
    
    % Patches
    for k=1:n
        patches.patches{k}.timeSolver = N;
        patches.patches{k}.timeOrder = 2;
    end
    
    %% Mass and stifness matrices and sollicitation vectors
    % Global
    glob.M = calc_mass(glob.S);
    glob.A = calc_rigi(glob.S);
    for k=1:n
        glob.M_in{k} = calc_mass(glob.S,'selgroup',getnumgroupelemwithparam(glob.S,'partition',k));
        glob.A_in{k} = calc_rigi(glob.S,'selgroup',getnumgroupelemwithparam(glob.S,'partition',k));
    end
    b_out = surfload(keepgroupelem(glob.S,getnumgroupelemwithparam(glob.S,'partition',0)),L2,'FX',-1);
    glob.b_out = b_out*loadFunction(glob.timeSolver);
    
    % Complementary subdomain
    globOut.M = calc_mass(globOut.S);
    globOut.A = calc_rigi(globOut.S);
    b = surfload(globOut.S,L2,'FX',-1);
    globOut.b = b*loadFunction(globOut.timeSolver);
    
    % Patches
    for k=1:n
        patches.patches{k}.M = calc_mass(patches.patches{k}.S);
        patches.patches{k}.A = calc_rigi(patches.patches{k}.S);
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
    
    %% Save variables
    save(fullfile(pathname,'problem.mat'),'glob','globOut','patches','interfaces','N','D','D_patch','L1','L2');
else
    load(fullfile(pathname,'problem.mat'),'glob','globOut','patches','interfaces','N','D','D_patch','L1','L2');
end

%% Direct solver
if directSolver
    DS = DirectSolver();
    DS.changeOfVariable = false;
    DS.timeSolver = N;
    DS.timeOrder = 2;
    DS.display = true;
    
    [Ut_ref,wt_ref,lambdat_ref,outputt_ref,vUt_ref,vwt_ref,vlambdat_ref,aUt_ref,awt_ref,alambdat_ref] = DS.solve(globOut,patches,interfaces);
    vUt_ref = unfreevector(globOut.S,vUt_ref)-calc_init_dirichlet(globOut.S);
    aUt_ref = unfreevector(globOut.S,aUt_ref)-calc_init_dirichlet(globOut.S);
    for k=1:n
        vwt_ref{k} = unfreevector(patches.patches{k}.S,vwt_ref{k})-calc_init_dirichlet(patches.patches{k}.S);
        vlambdat_ref{k} = unfreevector(interfaces.interfaces{k}.S,vlambdat_ref{k})-calc_init_dirichlet(interfaces.interfaces{k}.S);
        awt_ref{k} = unfreevector(patches.patches{k}.S,awt_ref{k})-calc_init_dirichlet(patches.patches{k}.S);
        alambdat_ref{k} = unfreevector(interfaces.interfaces{k}.S,alambdat_ref{k})-calc_init_dirichlet(interfaces.interfaces{k}.S);
    end
    save(fullfile(pathname,'reference_solution.mat'),'Ut_ref','wt_ref','lambdat_ref','outputt_ref','vUt_ref','vwt_ref','vlambdat_ref','aUt_ref','awt_ref','alambdat_ref');
else
    load(fullfile(pathname,'reference_solution.mat'),'Ut_ref','wt_ref','lambdat_ref','outputt_ref','vUt_ref','vwt_ref','vlambdat_ref','aUt_ref','awt_ref','alambdat_ref');
end

%% Outputs
fprintf('\n')
fprintf('spatial dimension = %d for U_ref\n',size(Ut_ref,1))
for k=1:n
    fprintf('                  = %d for w_ref{%u}\n',size(wt_ref{k},1),k)
    fprintf('                  = %d for lambda_ref{%u}\n',size(lambdat_ref{k},1),k)
end
fprintf('time solver : %s\n',class(N));
fprintf('nb time steps = %g\n',getnt(N))
fprintf('nb time dofs  = %g\n',getnbtimedof(N))
fprintf('elapsed time = %f s\n',outputt_ref.time)

%% Global-local Iterative solver
if iterativeSolver
    IS = IterativeSolver();
    IS.maxIterations = 50;
    IS.tolerance = eps;
    IS.relaxation = 'Aitken';
    IS.updateRelaxationParameter = true;
    IS.errorCriterion = 'reference';
    IS.referenceSolution = {Ut_ref,wt_ref,lambdat_ref};
    IS.display = true;
    IS.displayIterations = true;
    
    [Ut,wt,lambdat,outputt,vUt,vwt,vlambdat,aUt,awt,alambdat] = IS.solve(glob,patches,interfaces);
    vUt = unfreevector(glob.S,vUt)-calc_init_dirichlet(glob.S);
    aUt = unfreevector(glob.S,aUt)-calc_init_dirichlet(glob.S);
    for k=1:n
        vwt{k} = unfreevector(patches.patches{k}.S,vwt{k})-calc_init_dirichlet(patches.patches{k}.S);
        vlambdat{k} = unfreevector(interfaces.interfaces{k}.S,vlambdat{k})-calc_init_dirichlet(interfaces.interfaces{k}.S);
        awt{k} = unfreevector(patches.patches{k}.S,awt{k})-calc_init_dirichlet(patches.patches{k}.S);
        alambdat{k} = unfreevector(interfaces.interfaces{k}.S,alambdat{k})-calc_init_dirichlet(interfaces.interfaces{k}.S);
    end
    save(fullfile(pathname,'solution.mat'),'Ut','wt','lambdat','outputt','vUt','vwt','vlambdat','aUt','awt','alambdat');
else
    load(fullfile(pathname,'solution.mat'),'Ut','wt','lambdat','outputt','vUt','vwt','vlambdat','aUt','awt','alambdat');
end

%% Outputs
fprintf('\n')
fprintf('spatial dimension = %d for U\n',size(Ut,1))
for k=1:n
    fprintf('                  = %d for w{%u}\n',size(wt{k},1),k)
    fprintf('                  = %d for lambda{%u}\n',size(lambdat{k},1),k)
end
fprintf('nb time steps = %g\n',getnt(N))
fprintf('nb time dofs  = %g\n',getnbtimedof(N))
fprintf('elapsed time = %f s\n',outputt.totalTime)

%% Display
if displaySolution
    %% Display domains and meshes
    % plotDomain(D,D_patch);
    % mysaveas(pathname,'domain_global_patches',formats,renderer);
    % mymatlab2tikz(pathname,'domain_global_patches.tex');
    
    numbers = getnumber(patches);
    figure('Name',['Complement subdomain and patches #' num2str([numbers{:}])])
    clf
    h1 = plot(D,'FaceColor',getfacecolor(1));
    hold on
    h2 = plot(L1,'EdgeColor',getfacecolor(5));
    h3 = plot(L2,'EdgeColor',getfacecolor(6));
    h_patch = cell(1,n);
    leg_patch = cell(1,n);
    for k=1:n
        h_patch{k} = plot(D_patch{k},'FaceColor',getfacecolor(1+k));
        leg_patch{k} = ['$\Lambda_{' num2str(k) '}$'];
    end
    hold off
    set(gca,'FontSize',fontsize)
    l = legend([h1(1),h2(1),h3(1),h_patch{:}],'$\Omega \setminus \Lambda$','$\Gamma_D$','$\Gamma_N$',leg_patch{:},'Location','NorthEastOutside');
    set(l,'Interpreter','latex')
    axis image
    axis off
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
    h2 = plot(L1,'EdgeColor',getfacecolor(5));
    h3 = plot(L2,'EdgeColor',getfacecolor(6));
    h_patch = cell(1,n);
    leg_patch = cell(1,n);
    for k=1:n
        h_patch{k} = plot(patches.patches{k}.S,'FaceColor',getfacecolor(1+k));
        leg_patch{k} = ['$\Lambda_{' num2str(k) '}$'];
    end
    hold off
    set(gca,'FontSize',fontsize)
    % l = legend([h1(1),h2(1),h3(1),h_patch{:}],...
    %     '$\Omega \setminus \Lambda$','$\Gamma_D$','$\Gamma_N$',leg_patch{:},'Location','NorthEastOutside');
    % set(l,'Interpreter','latex')
    mysaveas(pathname,'mesh_global_patches',formats,renderer);
    
    %% Display evolutions of error indicator, stagnation indicator, CPU time w.r.t. number of iterations
    plotError(outputt);
    mysaveas(pathname,'error','fig');
    mymatlab2tikz(pathname,'error.tex');
    
    plotStagnation(outputt);
    mysaveas(pathname,'stagnation','fig');
    mymatlab2tikz(pathname,'stagnation.tex');
    
    plotErrorGlobalSolution(outputt);
    mysaveas(pathname,'error_global_solution','fig');
    mymatlab2tikz(pathname,'error_global_solution.tex');
    
    plotStagnationGlobalSolution(outputt);
    mysaveas(pathname,'stagnation_global_solution','fig');
    mymatlab2tikz(pathname,'stagnation_global_solution.tex');
    
    plotCPUTime(outputt,'legend',false);
    mysaveas(pathname,'cpu_time','fig');
    mymatlab2tikz(pathname,'cpu_time.tex');
    
    plotRelaxationParameter(outputt,'legend',false);
    mysaveas(pathname,'relaxation_parameter','fig');
    mymatlab2tikz(pathname,'relaxation_parameter.tex');
    
    %% Display evolutions of solutions
    i = 1;
    % for i=1:2
        % evolAllSolutions(glob,patches,interfaces,Ut,wt,lambdat,'displ',i,'filename',['all_solutions_' num2str(i)]],'pathname',pathname);
        % evolAllSolutions(glob,patches,interfaces,Ut,wt,lambdat,'displ',i,'surface',true,'filename',['all_solutions_' num2str(i) '_surface'],'pathname',pathname);
        % evolAllSolutions(glob,patches,interfaces,vUt,vwt,vlambdat,'displ',i,'filename',['all_velocities_' num2str(i)]],'pathname',pathname);
        % evolAllSolutions(glob,patches,interfaces,vUt,vwt,vlambdat,'displ',i,'surface',true,'filename',['all_velocities_' num2str(i) '_surface'],'pathname',pathname);
        % evolAllSolutions(glob,patches,interfaces,aUt,awt,alambdat,'displ',i,'filename',['all_accelerations_' num2str(i)]],'pathname',pathname);
        % evolAllSolutions(glob,patches,interfaces,aUt,awt,alambdat,'displ',i,'surface',true,'filename',['all_accelerations_' num2str(i) '_surface'],'pathname',pathname);
        
        evolGlobalSolution(glob,Ut,'displ',i,'filename',['global_solution_' num2str(i)],'pathname',pathname);
        evolGlobalSolution(glob,vUt,'displ',i,'filename',['global_solution_velocity_' num2str(i)],'pathname',pathname);
        evolGlobalSolution(glob,aUt,'displ',i,'filename',['global_solution_acceleration_' num2str(i)],'pathname',pathname);
        
        % evolLocalSolution(patches,wt,'displ',i,'filename',['local_solution_' num2str(i)],'pathname',pathname);
        % evolLocalSolution(patches,vwt,'displ',i,'filename',['local_solution_velocity_' num2str(i)],'pathname',pathname);
        % evolLocalSolution(patches,awt,'displ',i,'filename',['local_solution_acceleration_' num2str(i)],'pathname',pathname);
        
        % evolLagrangeMultiplier(interfaces,lambdat,'displ',i,'filename',['Lagrange_multiplier_' num2str(i)],'pathname',pathname);
        % evolLagrangeMultiplier(interfaces,vlambdat,'displ',i,'filename',['Lagrange_multiplier_velocity_' num2str(i)],'pathname',pathname);
        % evolLagrangeMultiplier(interfaces,alambdat,'displ',i,'filename',['Lagrange_multiplier_acceleration_' num2str(i)],'pathname',pathname);
        
        evolMultiscaleSolution(glob,patches,interfaces,Ut,wt,'displ',i,'filename',['multiscale_solution_' num2str(i)],'pathname',pathname);
        evolMultiscaleSolution(glob,patches,interfaces,vUt,vwt,'displ',i,'filename',['multiscale_solution_velocity_' num2str(i)],'pathname',pathname);
        evolMultiscaleSolution(glob,patches,interfaces,aUt,awt,'displ',i,'filename',['multiscale_solution_acceleration_' num2str(i)],'pathname',pathname);
        
        evolGlobalLocalSolution(glob,patches,interfaces,Ut,wt,'displ',i,'filename',['global_local_solution_' num2str(i)],'pathname',pathname);
        evolGlobalLocalSolution(glob,patches,interfaces,vUt,vwt,'displ',i,'filename',['global_local_solution_velocity_' num2str(i)],'pathname',pathname);
        evolGlobalLocalSolution(glob,patches,interfaces,aUt,awt,'displ',i,'filename',['global_local_solution_acceleration_' num2str(i)],'pathname',pathname);
    % end
    
    %% Display solutions at different instants
    [t,rep] = gettevol(N);
    for k=1:floor(length(rep)/5):length(rep)
        close all
        Uk = getmatrixatstep(Ut,rep(k));
        wk = cellfun(@(x) getmatrixatstep(x,rep(k)),wt,'UniformOutput',false);
        lambdak = cellfun(@(x) getmatrixatstep(x,rep(k)),lambdat,'UniformOutput',false);
        vUk = getmatrixatstep(vUt,rep(k));
        vwk = cellfun(@(x) getmatrixatstep(x,rep(k)),vwt,'UniformOutput',false);
        vlambdak = cellfun(@(x) getmatrixatstep(x,rep(k)),vlambdat,'UniformOutput',false);
        aUk = getmatrixatstep(aUt,rep(k));
        awk = cellfun(@(x) getmatrixatstep(x,rep(k)),awt,'UniformOutput',false);
        alambdak = cellfun(@(x) getmatrixatstep(x,rep(k)),alambdat,'UniformOutput',false);
        
        i = 1;
        % for i=1:2
            % plotAllSolutions(glob,patches,interfaces,Uk,wk,lambdak,'displ',i);
            % mysaveas(pathname,['all_solutions_' num2str(i) '_t' num2str(k-1)],formats,renderer);
            % plotAllSolutions(glob,patches,interfaces,vUk,vwk,vlambdak,'displ',i);
            % mysaveas(pathname,['all_velocities_' num2str(i) '_t' num2str(k-1)],formats,renderer);
            % plotAllSolutions(glob,patches,interfaces,aUk,awk,alambdak,'displ',i);
            % mysaveas(pathname,['all_accelerations_' num2str(i) '_t' num2str(k-1)],formats,renderer);
            
            plotGlobalSolution(glob,Uk,'displ',i);
            mysaveas(pathname,['global_solution_' num2str(i) '_t' num2str(k-1)],formats,renderer);
            plotGlobalSolution(glob,vUk,'displ',i);
            mysaveas(pathname,['global_solution_velocity_' num2str(i) '_t' num2str(k-1)],formats,renderer);
            plotGlobalSolution(glob,aUk,'displ',i);
            mysaveas(pathname,['global_solution_acceleration_' num2str(i) '_t' num2str(k-1)],formats,renderer);
            
            % plotLocalSolution(patches,wk,'displ',i);
            % mysaveas(pathname,['local_solution_' num2str(i) '_t' num2str(k-1)],formats,renderer);
            % plotLocalSolution(patches,vwk,'displ',i);
            % mysaveas(pathname,['local_solution_velocity_' num2str(i) '_t' num2str(k-1)],formats,renderer);
            % plotLocalSolution(patches,awk,'displ',i);
            % mysaveas(pathname,['local_solution_acceleration_' num2str(i) '_t' num2str(k-1)],formats,renderer);
            
            % plotLagrangeMultiplier(interfaces,lambdak,'displ',i);
            % mysaveas(pathname,['Lagrange_multiplier_' num2str(i) '_t' num2str(k-1)],formats,renderer);
            % plotLagrangeMultiplier(interfaces,vlambdak,'displ',i);
            % mysaveas(pathname,['Lagrange_multiplier_velocity_' num2str(i) '_t' num2str(k-1)],formats,renderer);
            % plotLagrangeMultiplier(interfaces,alambdak,'displ',i);
            % mysaveas(pathname,['Lagrange_multiplier_acceleration_' num2str(i) '_t' num2str(k-1)],formats,renderer);
            
            plotMultiscaleSolution(glob,patches,interfaces,Uk,wk,'displ',i);
            mysaveas(pathname,['multiscale_solution_' num2str(i) '_t' num2str(k-1)],formats,renderer);
            plotMultiscaleSolution(glob,patches,interfaces,vUk,vwk,'displ',i);
            mysaveas(pathname,['multiscale_solution_velocity_' num2str(i) '_t' num2str(k-1)],formats,renderer);
            plotMultiscaleSolution(glob,patches,interfaces,aUk,awk,'displ',i);
            mysaveas(pathname,['multiscale_solution_acceleration_' num2str(i) '_t' num2str(k-1)],formats,renderer);
            
            plotGlobalLocalSolution(glob,patches,interfaces,Uk,wk,'displ',i);
            mysaveas(pathname,['global_local_solution_' num2str(i) '_t' num2str(k-1)],formats,renderer);
            plotGlobalLocalSolution(glob,patches,interfaces,vUk,vwk,'displ',i);
            mysaveas(pathname,['global_local_solution_velocity_' num2str(i) '_t' num2str(k-1)],formats,renderer);
            plotGlobalLocalSolution(glob,patches,interfaces,aUk,awk,'displ',i);
            mysaveas(pathname,['global_local_solution_acceleration_' num2str(i) '_t' num2str(k-1)],formats,renderer);
        % end
    end
    
end

% myparallel('stop');
