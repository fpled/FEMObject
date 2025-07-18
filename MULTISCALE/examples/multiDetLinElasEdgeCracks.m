%% Multiscale deterministic linear elasticity problem with n edge cracks %%
%%-----------------------------------------------------------------------%%

% clc
clearvars
close all
% myparallel('start');

%% Input data
setProblem = true;
directSolver = true;
iterativeSolver = true;
displaySolution = true;

n = 1; % number of cracks n = 1
loading = 'Tension'; % 'Tension' or 'Shear'
filename = ['linElas' num2str(n) 'EdgeCracks' loading];
pathname = fullfile(getfemobjectoptions('path'),'MULTISCALE',...
    'examples','multiscaleDet',filename);
if ~exist(pathname,'dir')
    mkdir(pathname);
end
formats = {'fig','epsc'};
renderer = 'OpenGL';

%% Problem
if setProblem
    %% Domains and meshes
    % Global
    glob = Global();
    globOut = GlobalOutside();
    
    L = 16;
    w = 7;
    D = DOMAIN(2,[0.0,-L/2],[w,L/2]);
    
    option = 'DEFO'; % plane strain
    nbelem = [w*2,L*2];
    glob.S = build_model(D,'nbelem',nbelem,'option',option);
    % cl = 0.25;
    % glob.S = build_model(D,'cl',cl,'option',option,'filename',fullfile(pathname,'gmsh_domain'));
    
    % Patches
    patches = Patches(n);
    
    a = w/2;
    l = 1;
    D_patch = cell(1,n);
    B_patch = cell(1,n);
    P_patch = cell(1,n);
    switch n
        case 1
            P_patch{1} = [a,0.0];
            B_patch{1} = LIGNE([0.0,0.0],P_patch{1});
            D_patch{1} = DOMAIN(2,[0.0,-l],P_patch{1}+[l,l]);
        otherwise
            error('Wrong number of patches')
    end
    
    cl_patch_D = 0.25;
    cl_patch_B = 0.05;
    for k=1:n
        patches.patches{k}.S = gmshDomainWithSingleEdgeCrack(D_patch{k},B_patch{k},cl_patch_D,cl_patch_B,fullfile(pathname,['gmsh_patch_' num2str(k) '_single_edge_crack']));
        patches.patches{k}.S = setoption(patches.patches{k}.S,option);
    end
    
    % Partition of global mesh
    glob = partition(glob,D_patch);
    
    %% Materials
    % Poisson ratio
    NU = 0.3;
    % Thickness
    DIM3 = 1;
    % Density
    RHO = 1;
    % Young modulus
    E_out = 1;
    E_patch = cell(1,n);
    E_in = cell(1,n);
    for k=1:n
        patch = patches.patches{k};
        % E_patch(x) = 1 + f(x)
        % E_in(x)    = 1
        % with f(x) = 1 if ||x-c||_Inf < L
        %           = 0 if ||x-c||_Inf >= L
        % L = norm(getsize(D_patch{k}),Inf)/4;
        % c = getcenter(D_patch{k});
        % f = @(x) distance(x,c,Inf)<L;
        % E_patch{k} = FENODEFIELD(ones(patch.S.nbnode,1) + double(squeeze(f(patch.S.node))));
        E_patch{k} = 1;
        E_in{k} = 1;
    end
    
    % Complementary subdomain
    mat_out = ELAS_ISOT('E',E_out,'NU',NU,'RHO',RHO,'DIM3',DIM3);
    mat_out = setnumber(mat_out,0);
    glob.S = setmaterial(glob.S,mat_out,getnumgroupelemwithparam(glob.S,'partition',0));
    
    % Patches
    mat_patch = MATERIALS();
    for k=1:n
        mat_patch{k} = ELAS_ISOT('E',E_patch{k},'NU',NU,'RHO',RHO,'DIM3',DIM3);
        mat_patch{k} = setnumber(mat_patch{k},k);
        patches.patches{k}.S = setmaterial(patches.patches{k}.S,mat_patch{k});
    end
    
    % Fictitious patches
    mat_in = MATERIALS();
    for k=1:n
        mat_in{k} = ELAS_ISOT('E',E_in{k},'NU',NU,'RHO',RHO,'DIM3',DIM3);
        mat_in{k} = setnumber(mat_in{k},k);
        glob.S = setmaterial(glob.S,mat_in{k},getnumgroupelemwithparam(glob.S,'partition',k));
    end
    
    %% Dirichlet boundary conditions
    LU = LIGNE([0.0,L/2],[w,L/2]);
    LL = LIGNE([0.0,-L/2],[w,-L/2]);
    
    % Global
    glob.S = final(glob.S);
    switch lower(loading)
        case 'tension'
            glob.S = addcl(glob.S,POINT([a,0.0]),'UY');
            glob.S = addcl(glob.S,POINT([w,0.0]),'U');
        case 'shear'
            glob.S = addcl(glob.S,LL);
        otherwise
            error('Wrong loading case')
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
        patches.patches{k}.S = final(patches.patches{k}.S,'duplicate');
    end
    
    % Interfaces
    interfaces = Interfaces(patches,glob);
    
    %% Stiffness matrices and sollicitation vectors
    % Traction force density
    f = 1;
    
    % Global
    glob.A = calc_rigi(glob.S);
    for k=1:n
        glob.A_in{k} = calc_rigi(glob.S,'selgroup',getnumgroupelemwithparam(glob.S,'partition',k));
    end
    switch lower(loading)
        case 'tension'
            glob.b_out = surfload(keepgroupelem(glob.S,getnumgroupelemwithparam(glob.S,'partition',0)),LU,{'FX','FY'},[0;f]);
            glob.b_out = glob.b_out + surfload(keepgroupelem(glob.S,getnumgroupelemwithparam(glob.S,'partition',0)),LL,{'FX','FY'},[0;-f]);
        case 'shear'
            glob.b_out = surfload(keepgroupelem(glob.S,getnumgroupelemwithparam(glob.S,'partition',0)),LU,{'FX','FY'},[f;0]);
        otherwise
            error('Wrong loading case')
    end
    
    % Complementary subdomain
    globOut.A = calc_rigi(globOut.S);
    switch lower(loading)
        case 'tension'
            globOut.b = surfload(globOut.S,LU,{'FX','FY'},[0;f]);
            globOut.b = globOut.b + surfload(globOut.S,LL,{'FX','FY'},[0;-f]);
        case 'shear'
            globOut.b = surfload(globOut.S,LU,{'FX','FY'},[f;0]);
        otherwise
            error('Wrong loading case')
    end
    
    % Patches
    for k=1:n
        patches.patches{k}.A = calc_rigi(patches.patches{k}.S);
        patches.patches{k}.b = sparse(getnbddlfree(patches.patches{k}.S),1);
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
    glob.increment = true;
    
    % Local problems
    for k=1:n
        patches.patches{k}.changeOfVariable = false;
        patches.patches{k}.increment = true;
    end
    
    %% Save variables
    save(fullfile(pathname,'problem.mat'),'glob','globOut','patches','interfaces','D','D_patch','B_patch');
else
    load(fullfile(pathname,'problem.mat'),'glob','globOut','patches','interfaces','D','D_patch','B_patch');
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
    IS.maxIterations = 100;
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
    plotDomain(D,cellfun(@(x,y) {x,y},D_patch,B_patch,'UniformOutput',false));
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
    ampl = 0;
    % ampl = getsize(glob.S)/max(abs(U))/20;
    
    for i=1:2
        % plotAllSolutions(glob,patches,interfaces,U,w,lambda,'displ',i,'ampl',ampl);
        % mysaveas(pathname,['all_solutions_' num2str(i)],formats,renderer);
        
        plotGlobalSolution(glob,U,'displ',i,'ampl',ampl);
        mysaveas(pathname,['global_solution_' num2str(i)],formats,renderer);
        
        % plotLocalSolution(patches,w,'displ',i,'ampl',ampl);
        % mysaveas(pathname,['local_solution_' num2str(i)],formats,renderer);
        
        % plotLagrangeMultiplier(interfaces,lambda,'displ',i,'ampl',ampl);
        % mysaveas(pathname,['Lagrange_multiplier_' num2str(i)],formats,renderer);
        
        plotMultiscaleSolution(glob,patches,interfaces,U,w,'displ',i,'ampl',ampl);
        mysaveas(pathname,['multiscale_solution_' num2str(i)],formats,renderer);
        
        plotGlobalLocalSolution(glob,patches,interfaces,U,w,'displ',i,'ampl',ampl);
        mysaveas(pathname,['global_local_solution_' num2str(i)],formats,renderer);
    end
end

% myparallel('stop');
