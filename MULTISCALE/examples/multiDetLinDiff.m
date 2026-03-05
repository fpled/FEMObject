%% Multiscale deterministic linear diffusion problem %%
%%---------------------------------------------------%%

% clc
clearvars
close all
% myparallel('start');

%% Input data
setProblem = true;
directSolver = true;
iterativeSolver = true;
displaySolution = true;

Dim = 2; % space dimension Dim = 2, 3
n = 4; % number of patches n = 1, 2, 4
filename = ['linDiff' num2str(n) 'Patches_' num2str(Dim) 'D'];
pathname = fullfile(getfemobjectoptions('path'),'MULTISCALE',...
    'examples','multiscaleDet',filename);
if ~exist(pathname,'dir')
    mkdir(pathname);
end

formats = {'epsc'};
renderer = 'OpenGL';

%% Problem
if setProblem
    %% Domains and meshes
    % Global
    glob = Global();
    globOut = GlobalOutside();
    
    if Dim==2
        D = DOMAIN(2,[0.0,0.0],[1.0,1.0]);
    elseif Dim==3
        D = DOMAIN(3,[0.0,0.0,0.0],[1.0,1.0,1.0]);
    end
    
    nbelem = repmat(20,1,Dim);
    glob.S = build_model(D,'nbelem',nbelem);
    % cl = 0.05;
    % glob.S = build_model(D,'cl',cl,'filename',fullfile(pathname,'gmsh_domain'));
    
    % Patches
    patches = Patches(n);
    
    D_patch = cell(1,n);
    switch n
        case 1
            if Dim==2
                D_patch{1} = DOMAIN(2,[0.4,0.4],[0.6,0.6]);
            elseif Dim==3
                D_patch{1} = DOMAIN(3,[0.4,0.4,0.4],[0.6,0.6,0.6]);
            end
        case 2
            if Dim==2
                D_patch{1} = DOMAIN(2,[0.1,0.1],[0.3,0.3]);
                D_patch{2} = DOMAIN(2,[0.7,0.7],[0.9,0.9]);
            elseif Dim==3
                D_patch{1} = DOMAIN(3,[0.1,0.1,0.1],[0.3,0.3,0.3]);
                D_patch{2} = DOMAIN(3,[0.7,0.7,0.7],[0.9,0.9,0.9]);
            end
        case 4
            if Dim==2
                D_patch{1} = DOMAIN(2,[0.1,0.1],[0.3,0.3]);
                D_patch{2} = DOMAIN(2,[0.1,0.7],[0.3,0.9]);
                D_patch{3} = DOMAIN(2,[0.7,0.7],[0.9,0.9]);
                D_patch{4} = DOMAIN(2,[0.7,0.1],[0.9,0.3]);
            elseif Dim==3
                D_patch{1} = DOMAIN(3,[0.1,0.1,0.1],[0.3,0.3,0.3]);
                D_patch{2} = DOMAIN(3,[0.1,0.7,0.1],[0.3,0.9,0.3]);
                D_patch{3} = DOMAIN(3,[0.7,0.7,0.7],[0.9,0.9,0.9]);
                D_patch{4} = DOMAIN(3,[0.7,0.1,0.7],[0.9,0.3,0.9]);
            end
        otherwise
            error('Wrong number of patches')
    end
    
    nbelem_patch = repmat(40,Dim);
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
    % Linear diffusion coefficient
    K_out = 1;
    K_patch = cell(1,n);
    K_in = cell(1,n);
    for k=1:n
        patch = patches.patches{k};
        % K_patch(x) = 1 + beta_patch * f(x)
        % K_in(x)    = 1 + beta_in * f(x)
        % with f(x) = alpha*exp( -ampl*||x-c||_2^2/L^2 ) if ||x-c||_Inf < L
        %           = 0                                 if ||x-c||_Inf >= L
        % alpha = 10;
        % ampl = 2;
        % L = norm(getsize(D_patch{k}),Inf)/4;
        % c = getcenter(D_patch{k});
        % f = @(x) (distance(x,c,Inf)<L) * alpha * exp(-ampl*distance(x,c,2).^2/L^2);
        % beta_patch = 1;
        % K_patch{k} = FENODEFIELD(1 + beta_patch * squeeze(f(patch.S.node)));
        % beta_in = 0;
        % K_in{k} = FENODEFIELD(1 + beta_in * squeeze(f(glob.S.node)));
        
        % K_patch(x) = 1 + f(x)
        % K_in(x)    = 1
        % with f(x) = 1 if ||x-c||_Inf < L
        %           = 0 if ||x-c||_Inf >= L
        L = norm(getsize(D_patch{k}),Inf)/4;
        c = getcenter(D_patch{k});
        f = @(x) distance(x,c,Inf)<L;
        K_patch{k} = FENODEFIELD(ones(patch.S.nbnode,1) + double(squeeze(f(patch.S.node))));
        K_in{k} = 1;
    end
    
    % Complementary subdomain
    mat_out = FOUR_ISOT('k',K_out);
    mat_out = setnumber(mat_out,0);
    glob.S = setmaterial(glob.S,mat_out,getnumgroupelemwithparam(glob.S,'partition',0));
    
    % Patches
    mat_patch = MATERIALS();
    for k=1:n
        mat_patch{k} = FOUR_ISOT('k',K_patch{k});
        mat_patch{k} = setnumber(mat_patch{k},k);
        patches.patches{k}.S = setmaterial(patches.patches{k}.S,mat_patch{k});
    end
    
    % Fictitious patches
    mat_in = MATERIALS();
    for k=1:n
        mat_in{k} = FOUR_ISOT('k',K_in{k});
        mat_in{k} = setnumber(mat_in{k},k);
        glob.S = setmaterial(glob.S,mat_in{k},getnumgroupelemwithparam(glob.S,'partition',k));
    end
    
    %% Dirichlet boundary conditions
    % Global
    glob.S = final(glob.S);
    glob.S = addcl(glob.S,[]);
    glob.S_out = getfinalmodelpart(glob.S,0);
    % glob.S_in = cell(1,n);
    % for k=1:n
    %     glob.S_in{k} = getfinalmodelpart(glob.S,k);
    % end
    
    % Complementary subdomain
    globOut.S = glob.S_out;
    
    % Patches
    for k=1:n
        patches.patches{k}.S = final(patches.patches{k}.S);
    end
    
    % Interfaces
    interfaces = Interfaces(patches);
    
    %% Stiffness matrices and sollicitation vectors
    % Source term
    f = 100;
    
    % Global
    glob.A = calc_rigi(glob.S);
    for k=1:n
        glob.A_in{k} = calc_rigi(glob.S,'selgroup',getnumgroupelemwithparam(glob.S,'partition',k));
    end
    glob.b_out = bodyload(keepgroupelem(glob.S,getnumgroupelemwithparam(glob.S,'partition',0)),[],'QN',f);
    
    % Complementary subdomain
    globOut.A = calc_rigi(globOut.S);
    globOut.b = bodyload(globOut.S,[],'QN',f);
    
    % Patches
    for k=1:n
        patches.patches{k}.A = calc_rigi(patches.patches{k}.S);
        patches.patches{k}.b = bodyload(patches.patches{k}.S,[],'QN',f);
    end
    
    %% Mass matrices
    for k=1:n
        interfaces.interfaces{k}.M = calc_massgeom(interfaces.interfaces{k}.S);
    end
    
    %% Projection operators
    glob.P_out = calcProjection(glob);
    for k=1:n
        [interfaces.interfaces{k}.P_glob,numnode] = calcProjection(glob,interfaces.interfaces{k});
        % [interfaces.interfaces{k}.P_glob,numnode] = calcProjection(glob,interfaces.interfaces{k},'all','full',true);
        interfaces.interfaces{k}.P_globOut = interfaces.interfaces{k}.P_glob*glob.P_out';
        interfaces.interfaces{k}.P_patch = calcProjection(patches.patches{k},interfaces.interfaces{k});
        % plotProjectionOperator(glob,patches.patches{k},numnode);
    end
    
    %% Parameters for global and local problems
    % Global problem
    glob.increment = true;
    
    % Local problems
    for k=1:n
        patches.patches{k}.changeOfVariable = false;
        patches.patches{k}.increment = true;
    end
    
    %% Save variables
    save(fullfile(pathname,'problem.mat'),'glob','globOut','patches','interfaces','D','D_patch');
else
    load(fullfile(pathname,'problem.mat'),'glob','globOut','patches','interfaces','D','D_patch');
end

%% Direct solver
if directSolver
    DS = DirectSolver();
    DS.changeOfVariable = false;
    DS.display = true;
    
    [U_ref,w_ref,lambda_ref,output_ref] = DS.solve(globOut,patches,interfaces);
    save(fullfile(pathname,'reference_solution.mat'),'U_ref','w_ref','lambda_ref','output_ref');
else
    load(fullfile(pathname,'reference_solution.mat'),'U_ref','w_ref','lambda_ref','output_ref');
end

%% Outputs
fprintf('\n')
fprintf('spatial dimension = %d for U_ref\n',length(U_ref))
for k=1:n
    fprintf('                  = %d for w_ref{%u}\n',length(w_ref{k}),k)
    fprintf('                  = %d for lambda_ref{%u}\n',length(lambda_ref{k}),k)
end
fprintf('elapsed time = %f s\n',output_ref.time)

%% Global-local Iterative solver
if iterativeSolver
    IS = IterativeSolver();
    IS.maxIterations = 50;
    IS.tolerance = eps;
    IS.relaxation = 'Aitken';
    IS.updateRelaxationParameter = true;
    IS.errorCriterion = 'reference';
    IS.referenceSolution = {U_ref,w_ref,lambda_ref};
    IS.display = true;
    IS.displayIterations = true;
    
    [U,w,lambda,output] = IS.solve(glob,patches,interfaces);
    save(fullfile(pathname,'solution.mat'),'U','w','lambda','output');
else
    load(fullfile(pathname,'solution.mat'),'U','w','lambda','output');
end

%% Outputs
fprintf('\n')
fprintf('spatial dimension = %d for U\n',length(U))
for k=1:n
    fprintf('                  = %d for w{%u}\n',length(w{k}),k)
    fprintf('                  = %d for lambda{%u}\n',length(lambda{k}),k)
end
fprintf('elapsed time = %f s\n',output.totalTime)

%% Display
if displaySolution
    %% Display domains and meshes
    plotDomain(D,D_patch);
    mysaveas(pathname,'domain_global_patches',formats,renderer);
    mymatlab2tikz(pathname,'domain_global_patches.tex');
    
    % plotPartition(glob,'legend',false);
    % mysaveas(pathname,'mesh_partition',formats,renderer);
    
    plotModel(glob,patches,'legend',false);
    mysaveas(pathname,'mesh_global_patches',formats,renderer);
    
    % plotModel(glob);
    % plotModel(patches);
    % plotModel(interfaces);
    
    %% Display evolutions of error indicator, stagnation indicator, CPU time w.r.t. number of iterations
    plotError(output);
    mysaveas(pathname,'error','fig');
    mymatlab2tikz(pathname,'error.tex');
    
    plotStagnation(output);
    mysaveas(pathname,'stagnation','fig');
    mymatlab2tikz(pathname,'stagnation.tex');
    
    plotErrorGlobalSolution(output);
    mysaveas(pathname,'error_global_solution','fig');
    mymatlab2tikz(pathname,'error_global_solution.tex');
    
    plotStagnationGlobalSolution(output);
    mysaveas(pathname,'stagnation_global_solution','fig');
    mymatlab2tikz(pathname,'stagnation_global_solution.tex');
    
    plotCPUTime(output,'legend',false);
    mysaveas(pathname,'cpu_time','fig');
    mymatlab2tikz(pathname,'cpu_time.tex');
    
    plotRelaxationParameter(output,'legend',false);
    mysaveas(pathname,'relaxation_parameter','fig');
    mymatlab2tikz(pathname,'relaxation_parameter.tex');
    
    %% Display solutions
    % plotAllSolutions(glob,patches,interfaces,U,w,lambda);
    % mysaveas(pathname,'all_solutions',formats,renderer);
    
    plotGlobalSolution(glob,U);
    mysaveas(pathname,'global_solution',formats,renderer);
    
    % plotLocalSolution(patches,w);
    % mysaveas(pathname,'local_solution',formats,renderer);
    
    % plotLagrangeMultiplier(interfaces,lambda);
    % mysaveas(pathname,'Lagrange_multiplier',formats,renderer);
    
    plotMultiscaleSolution(glob,patches,interfaces,U,w);
    mysaveas(pathname,'multiscale_solution',formats,renderer);
    
    plotGlobalLocalSolution(glob,patches,interfaces,U,w);
    mysaveas(pathname,'global_local_solution',formats,renderer);
end

% myparallel('stop');
