classdef Global
    %
    
    %
    %   A MATLAB class for representing the global problem defined over
    %   the fictititous domain (complementary subdomain and fictitious
    %   patches)
    %
    
    properties
        S
        S_out
        A
        A_in
        b_out
        Atang
        M
        M_in
        C
        C_in
        b0_out
        u0
        v0
        P_out
        increment
        solver
        initializationType
        timeSolver
        timeOrder
        display
    end
    
    methods
        
        % Constructor
        function glob = Global(varargin)
            % Class Global
            %
            % function glob = Global(varargin)
            % glob.S: model associated to fictitious domain
            % glob.S_out: model associated to complementary subdomain
            % glob.A: stiffness operator associated to fictitious domain
            % glob.A_in: stiffness matrix associated to fictititous patch
            % glob.b_out: sollicitation vector associated to complementary subdomain
            % glob.Atang: tangent stiffness operator associated to fictitious domain
            % glob.M: mass operator associated to fictitious domain
            % glob.M_in: mass matrix associated to fictititous patch
            % glob.C: damping operator associated to fictitious domain
            % glob.C_in: damping matrix associated to fictititous patch
            % glob.b0_out: initial sollicitation vector associated to complementary subdomain, [] by default
            % glob.u0: initial solution associated to fictititous domain, [] by default
            % glob.v0: initial velocity associated to fictititous domain, [] by default
            % glob.P_out: projection operator from ficitious domain to complementary subdomain
            % glob.increment: reformulation of global problem on increments (instead of current iterates) (true or false), false by default
            % glob.solver: solver for global problem ([], NEWTONSOLVER), [] by default
            % if glob.solver is an instance of class NEWTONSOLVER,
            % glob.initializationType: type of initialization for iterative resolution of global problem ('zero', 'one' or 'previter'), 'zero' by default
            % glob.timeSolver: time solver for global problem ([], DGTIMESOLVER, EULERSOLVER, NEWMARKSOLVER), [] by default
            % glob.timeOrder: time order for global problem, 0 by default
            % glob.display: display (true or false), false by default
            
            switch nargin
                case 0
                    glob.S = [];
                    glob.S_out = [];
                    glob.A = [];
                    glob.A_in = [];
                    glob.b_out = [];
                    glob.Atang = [];
                    glob.M = [];
                    glob.M_in = [];
                    glob.C = [];
                    glob.C_in = [];
                    glob.b0_out = [];
                    glob.u0 = [];
                    glob.v0 = [];
                    glob.P_out = [];
                case 1
                    if isa(varargin{1},'Global')
                        glob = varargin{1};
                    else
                        glob.S = varargin{1};
                    end
                otherwise
                    glob.S = varargin{1};
                    glob.S_out = varargin{2};
                    glob.A = varargin{3};
                    glob.A_in = varargin{4};
                    glob.b_out = varargin{5};
                    glob.Atang = varargin{6};
                    glob.M = varargin{7};
                    glob.M_in = varargin{8};
                    glob.C = varargin{9};
                    glob.C_in = varargin{10};
                    glob.b0_out = varargin{11};
                    glob.u0 = varargin{12};
                    glob.v0 = varargin{13};
                    glob.P_out = varargin{14};
            end
            
            expectedInitializationTypes = {'zero','one','previter'};
            
            p = ImprovedInputParser;
            addParameter(p,'increment',false,@islogical);
            addParameter(p,'solver',[]);
            % addParameter(p,'initializationType','zero',@ischar);
            addParameter(p,'initializationType','zero',...
                @(x) any(validatestring(x,expectedInitializationTypes)))
            addParameter(p,'timeSolver',[]);
            addParameter(p,'timeOrder',0);
            addParameter(p,'display',false,@islogical);
            parse(p,varargin{:});
            glob = passMatchedArgsToProperties(p,glob);
            
        end
        
        function glob = partition(glob,patches)
            % function glob = partition(glob,patches)
            % Partition of elements of global mesh glob.S
            % 'partition' = 0: group of elements of mesh glob.S outside all patches (not in a strict sense)
            % 'partition' = k: group of elements of mesh glob.S inside patch #k (in a strict sense)
            
            n = numel(patches);
            
            glob.S = setparamgroupelem(glob.S,'partition',0);
            for k=1:n
                if isa(patches,'Patches')
                    [~,~,numelem] = intersect(glob.S,patches.patches{k}.S,'strict',1);
                elseif iscell(patches) && (isa(patches{k},'GEOMOBJECT') || isa(patches{k},'MODEL'))
                    [~,~,numelem] = intersect(glob.S,patches{k},'strict',1);
                else
                    error('Wrong input arguments')
                end
                [glob.S,newgroupelem] = separateelemwithnum(glob.S,numelem);
                glob.S = setparamgroupelem(glob.S,'partition',k,newgroupelem);
            end
            
            % numelem = cell(1,n);
            % parfor k=1:n
            %     [~,~,numelem{k}] = intersect(glob.S,patches.patches{k}.S,'strict',1);
            % end
            % for k=1:n
            %     [glob.S,newgroupelem] = separateelemwithnum(glob.S,numelem{k});
            %     glob.S = setparamgroupelem(glob.S,'partition',k,newgroupelem);
            % end
        end
        
        function [V,output,vV,aV] = solve(glob,interfaces,lambda_old,U_old,V_old,vU_old,vV_old,aU_old,aV_old)
            % function [V,output] = solve(glob,interfaces,lambda_old,U_old,V_old)
            % function [V,output,vV] = solve(glob,interfaces,lambda_old,U_old,V_old,vU_old,vV_old)
            % function [V,output,vV,aV] = solve(glob,interfaces,lambda_old,U_old,V_old,vU_old,vV_old,aU_old,aV_old)
            % Solves global problem defined over fictitious domain
            %
            % Inputs:
            % glob: Global
            % interfaces: Interfaces
            % lambda_old: 1-by-n cell of m_l-by-p double containing previous Lagrange multiplier lambda_{k-1}
            % U_old: m_U-by-p double containing previous global iterate U_{k-1}
            % V_old: m_U-by-p double containing previous global solution V_{k-1}
            % For first- and second-order time-dependent problems,
            % vU_old: m_U-by-p double containing previous global iterate velocity vU_{k-1}
            % vV_old: m_U-by-p double containing previous global solution velocity vV_{k-1}
            % For second-order time-dependent problems,
            % aU_old: m_U-by-p double containing previous global iterate acceleration aU_{k-1}
            % aV_old: m_U-by-p double containing previous global solution acceleration aV_{k-1}
            %
            % Outputs:
            % V: m_U-by-p double containing current global solution V_{k}
            % output.time : 1-by-1 double containing CPU time t
            % For time-dependent problems,
            % output.result: structure containing outputs of time solver
            % For first- and second-order time-dependent problems,
            % vV: m_U-by-p double containing current global solution velocity vV_{k}
            % For second-order time-dependent problems,
            % aV: m_U-by-p double containing current global solution acceleration aV_{k}
            % where
            % k is the current iteration number
            % m_U is the dimension of the spatial approximation space of global solution U
            % m_l is the dimension of the spatial approximation space of Lagrange multiplier lambda
            % p is the dimension of the time approximation space of global solution U and Lagrange multiplier lambda
            
            t = tic;
            
            if nargin<5 || isempty(V_old)
                glob.increment = false;
            elseif ~isempty(glob.timeSolver)
                if (nargin<7 || isempty(vV_old)) && glob.timeOrder>=1
                    glob.increment = false;
                end
                if (nargin<9 || isempty(aV_old)) && glob.timeOrder>=2
                    glob.increment = false;
                end
            end
            
            glob = initializeRightHandSide(glob);
            if isempty(glob.timeSolver)
                b = setRightHandSide(glob,interfaces,lambda_old,U_old);
            else
                if glob.timeOrder==1
                    b = setRightHandSide(glob,interfaces,lambda_old,U_old,vU_old);
                 elseif glob.timeOrder==2
                    b = setRightHandSide(glob,interfaces,lambda_old,U_old,vU_old,aU_old);
                end
            end
            
            if glob.increment && (isempty(glob.solver) || ~isanlsolver(glob.solver))
                b = b - glob.A*V_old;
                if ~isempty(glob.timeSolver)
                    if glob.timeOrder==1
                        b = b - glob.M*vV_old;
                    elseif glob.timeOrder==2
                        b = b - glob.M*aV_old;
                        if ~isempty(glob.C)
                            b = b - glob.C*vV_old;
                        end
                    end
                end
            end
            
            % Global soluton V
            if isempty(glob.timeSolver)
                % Time-independent problem
                if isempty(glob.solver)
                    V = glob.A\b;
                else
                    if isa(glob.solver,'NEWTONSOLVER')
                        V0 = setInitialSolution(glob.initializationType,b,V_old);
                        V = solve(glob.solver,b,glob.A,glob.Atang,V0);
                    end
                end
            else
                % Time-dependent problem
                if glob.timeOrder==1
                    [V,output.result,vV] = dsolve(glob.timeSolver,b,glob.M,glob.A,glob.u0);
                elseif glob.timeOrder==2
                    [V,output.result,vV,aV] = ddsolve(glob.timeSolver,b,glob.M,glob.A,glob.C,glob.u0,glob.v0);
                end
            end
            
            if glob.increment && (isempty(glob.solver) || ~isanlsolver(glob.solver))
                V = V_old + V;
                if ~isempty(glob.timeSolver)
                    if glob.timeOrder>=1
                        vV = vV_old + vV;
                    end
                    if glob.timeOrder>=2
                        aV = aV_old + aV;
                    end
                end
            end
            
            output.time = toc(t);
            
            if glob.display
                fprintf('\nElapsed time = %f s\n',output.time);
            end
            
        end
        
        function glob = linearizeOperator(glob,U)
            % function glob = linearizeOperator(glob,U)
            % Linearize stiffness operator as the tangent stiffness operator evaluated at global solution U
            % glob: Global
            % U: global solution
            
            if ~isempty(glob.solver) && isanlsolver(glob.solver)
                glob.solver = [];
                glob.A = glob.Atang(U);
            end
            
        end
        
        function glob = initializeRightHandSide(glob)
            % function glob = initializeRightHandSide(glob)
            % Initialize right-hand side of global problem
            
            if isempty(glob.b_out)
                sz = getnbddlfree(glob.S);
                glob.b_out = sparse(sz,1);
                if ~isempty(glob.timeSolver)
                    T = gettimemodel(glob.timeSolver);
                    glob.b_out = glob.b_out*zero(T);
                end
            end
        end
        
        function b = setRightHandSide(glob,interfaces,lambda_old,U_old,vU_old,aU_old)
            % function b = setRightHandSide(glob,interfaces,lambda_old,U_old)
            % function b = setRightHandSide(glob,interfaces,lambda_old,U_old,vU_old)
            % function b = setRightHandSide(glob,interfaces,lambda_old,U_old,vU_old,aU_old)
            % Set right-hand side of global problem
            % glob: Global
            % interfaces: Interfaces
            % lambda_old: 1-by-n cell of m_l-by-p double containing previous Lagrange multiplier lambda_{k-1}
            % U_old: m_U-by-p double containing previous global iterate U_{k-1}
            % For first- and second-order time-dependent problems,
            % vU_old: m_U-by-p double containing previous global iterate velocity vU_{k-1}
            % For second-order time-dependent problems,
            % aU_old: m_U-by-p double containing previous global iterate acceleration aU_{k-1}
            
            b = glob.b_out;
            n = numel(interfaces);
            interface = interfaces.interfaces;
            for k=1:n
                B_glob = interface{k}.P_glob'*interface{k}.M;
                b = b - B_glob*lambda_old{k} + glob.A_in{k}*U_old;
                if ~isempty(glob.timeSolver)
                    if glob.timeOrder==1
                        b = b + glob.M_in{k}*vU_old;
                    elseif glob.timeOrder==2
                        b = b + glob.M_in{k}*aU_old;
                        if ~isempty(glob.C)
                            b = b + glob.C_in{k}*vU_old;
                        end
                    end
                end
            end
        end
        
        function [U0,V0] = setInitialCondition(glob)
            % function U0 = setInitialCondition(glob)
            % Sets global initial solution vector U0
            % function [U0,V0] = setInitialCondition(glob)
            % Sets global initial solution vector U0 and velocity vector V0
            % glob: Global
            
            U0 = initializeSolution(glob);
            if nargout>1
                V0 = initializeVelocity(glob);
            end
        end
        
        function U0 = initializeSolution(glob)
            % function U0 = initializeSolution(glob)
            
            sz = getnbddlfree(glob.S);
            if isempty(glob.u0)
                U0 = zeros(sz,1);
            else
                U0 = glob.u0;
            end
        end
        
        function V0 = initializeVelocity(glob)
            % function V0 = initializeVelocity(glob)
            
            sz = getnbddlfree(glob.S);
            if isempty(glob.v0)
                V0 = zeros(sz,1);
            else
                V0 = glob.v0;
            end
        end
        
        function [P,numnode] = calcProjection(glob,varargin)
            % function [P,numnode] = calcProjection(glob,varargin)
            % Calculates projection operator from model glob.S to model glob.S_out if nargin==1
            % Calculates projection operator from model glob.S to model varargin{1}.S if varargin{1} is an instance of type GlobalOutside or Interface
            % Calculates projection operator from model glob.S to model varargin{1} if varargin{1} is an instance of class MODEL
            % glob: Global
            
            switch nargin
                case 1
                    if isempty(glob.S_out)
                        S_out = getmodelpart(glob.S,0);
                    else
                        S_out = glob.S_out;
                    end
                    [P,numnode] = calc_P_free(glob.S,S_out);
                otherwise
                    if isa(varargin{1},'GlobalOutside')
                        S_out = varargin{1}.S;
                        free = getcharin('free',varargin,true);
                        varargin = delcharin('free',varargin);
                        if free
                            [P,numnode] = calc_P_free(glob.S,S_out,varargin{2:end});
                        else
                            [P,numnode] = calc_P(glob.S,S_out,varargin{2:end});
                        end
                    elseif isa(varargin{1},'Interface')
                        S_interface = varargin{1}.S;
                        [P,numnode] = calcProjection(S_interface,glob.S,varargin{2:end});
                        P = P';
                    elseif isa(varargin{1},'MODEL')
                        [P,numnode] = calcProjection(glob.S,varargin{:});
                    else
                        if isempty(glob.S_out)
                            S_out = getmodelpart(glob.S,0);
                        else
                            S_out = glob.S_out;
                        end
                        free = getcharin('free',varargin,true);
                        varargin = delcharin('free',varargin);
                        if free
                            [P,numnode] = calc_P_free(glob.S,S_out,varargin{:});
                        else
                            [P,numnode] = calc_P(glob.S,S_out,varargin{:});
                        end
                    end
            end
        end
        
        function plotModel(glob,varargin)
            % function plotModel(glob)
            % Display model glob.S
            %
            % function plotModel(glob,k)
            % Display selected parts of model glob.S
            %
            % function plotModel(glob,'out')
            % Display part of model glob.S outside all patches
            %
            % function plotModel(glob,'in')
            % Display parts of model glob.S inside all patches
            %
            % function plotModel(glob,patch)
            % Display model glob.S as well as model patch.S
            %
            % function plotModel(glob,patches)
            % Display model glob.S as well as models patches.S of all patches
            %
            % function plotModel(glob,interface)
            % Display model glob.S as well as model interface.S
            %
            % function plotModel(glob,interfaces)
            % Display model glob.S as well as models interfaces.S of all interfaces
            
            p = ImprovedInputParser;
            addParameter(p,'legend',true,@islogical);
            addParameter(p,'FontSize',16,@isscalar);
            addParameter(p,'LineWidth',0.5,@isscalar);
            addParameter(p,'Interpreter','latex',@ischar);
            if isempty(varargin) || (ischar(varargin{1}) && ~strcmpi(varargin{1},'out') && ~strcmpi(varargin{1},'in'))
                parse(p,varargin{:})
            else
                parse(p,varargin{2:end})
            end
            
            options = varargin;
            varargin = delcharin({'legend','FontSize','LineWidth','Interpreter'},varargin);
            nargin = length(varargin)+1;
            
            part = getelementparam(glob.S,'partition');
            nbparts = length(unique([part{:}]));
            
            switch nargin
                case 1
                    figure('Name','Mesh of fictitious domain')
                    % set(gcf,'Name','Mesh of fictitious domain')
                    clf
                    hg = cell(1,nbparts);
                    leg = cell(1,nbparts);
                    for k=1:nbparts
                        S_part = getmodelpart(glob.S,k-1);
                        h = plot(S_part,'FaceColor',getfacecolor(k),'LineWidth',p.Results.LineWidth,varargin{:});
                        hg{k} = hggroup;
                        set(h(:),'Parent',hg{k});
                        if k==1
                            leg{k} = '$\Omega \setminus \Lambda$';
                        else
                            leg{k} = ['$\widetilde{\Lambda}_{' num2str(k-1) '}$'];
                        end
                    end
                    if p.Results.legend
                        legend([hg{:}],leg{:},'Location','NorthEastOutside','Interpreter',p.Results.Interpreter)
                    end
                    set(gca,'FontSize',p.Results.FontSize)
                case 2
                    if isa(varargin{1},'double')
                        if varargin{1}==0
                            figure('Name','Mesh of complementary subdomain')
                            % set(gcf,'Name','Mesh of complementary subdomain')
                            clf
                            % S_out = getmodelpart(glob.S,0);
                            h = plot(glob.S_out,'FaceColor',getfacecolor(1),'LineWidth',p.Results.LineWidth,varargin{:});
                            hg = hggroup;
                            set(h(:),'Parent',hg);
                            if p.Results.legend
                                legend(hg,'$\Omega \setminus \Lambda$','Location','NorthEastOutside','Interpreter',p.Results.Interpreter)
                            end
                            set(gca,'FontSize',p.Results.FontSize)
                        elseif numel(varargin{1})==1
                            k = varargin{1};
                            figure('Name',['Mesh of fictitious patch #' num2str(k)])
                            % set(gcf,'Name',['Mesh of fictitious patch #' num2str(k)])
                            clf
                            S_in = getmodelpart(glob.S,k);
                            h = plot(S_in,'FaceColor',getfacecolor(k+1),'LineWidth',p.Results.LineWidth,varargin{:});
                            hg = hggroup;
                            set(h(:),'Parent',hg);
                            if p.Results.legend
                                legend(hg,['$\widetilde{\Lambda}_{' num2str(k) '}$'],...
                                    'Location','NorthEastOutside','Interpreter',p.Results.Interpreter)
                            end
                            set(gca,'FontSize',p.Results.FontSize)
                        else
                            parts = varargin{1};
                            if ismember(0,parts)
                                k = setdiff(parts,0);
                                figname = ['Meshes of complementary subdomain and fictitious patches #' num2str(k)];
                            else
                                k = parts;
                                figname = ['Meshes of selected fictitious patches #' num2str(k)];
                            end
                            figure('Name',figname)
                            % set(gcf,'Name',figname)
                            clf
                            hg = cell(1,numel(parts));
                            leg = cell(1,numel(parts));
                            for i=1:numel(parts)
                                k = parts(i);
                                S_part = getmodelpart(glob.S,k);
                                h = plot(S_part,'FaceColor',getfacecolor(k+1),'LineWidth',p.Results.LineWidth,varargin{:});
                                hg{k} = hggroup;
                                set(h(:),'Parent',hg{k});
                                if k==0
                                    leg{i} = '$\Omega \setminus \Lambda$';
                                else
                                    leg{i} = ['$\widetilde{\Lambda}_{' num2str(k) '}$'];
                                end
                            end
                            if p.Results.legend
                                legend([hg{:}],leg{:},'Location','NorthEastOutside','Interpreter',p.Results.Interpreter)
                            end
                            set(gca,'FontSize',p.Results.FontSize)
                        end
                    elseif ischar(varargin{1})
                        if strcmp(varargin{1},'out')
                            options = delonlycharin('out',options);
                            plotModel(glob,0,options{:});
                        elseif strcmp(varargin{1},'in')
                            options = delonlycharin('in',options);
                            plotModel(glob,1:nbparts-1,options{:});
                        end
                    elseif isa(varargin{1},'Patch')
                        patch = varargin{1};
                        
                        figure('Name',['Meshes of complementary subdomain and patch #' num2str(patch.number)])
                        % set(gcf,'Name',['Meshes of complementary subdomain and patch #' num2str(patch.number)])
                        clf
                        % S_out = getmodelpart(glob.S,0);
                        h_out = plot(glob.S_out,'FaceColor',getfacecolor(1),'LineWidth',p.Results.LineWidth,varargin{:});
                        h_patch = plot(patch.S,'FaceColor',getfacecolor(patch.number+1),'LineWidth',p.Results.LineWidth,varargin{:});
                        hg_out = hggroup;
                        hg_patch = hggroup;
                        set(h_out(:),'Parent',hg_out);
                        set(h_patch(:),'Parent',hg_patch);
                        if p.Results.legend
                            legend([hg_out,hg_patch],'$\Omega \setminus \Lambda$',['$\Lambda_{' num2str(patch.number) '}$'],...
                                'Location','NorthEastOutside','Interpreter',p.Results.Interpreter)
                        end
                        set(gca,'FontSize',p.Results.FontSize)
                    elseif isa(varargin{1},'Patches')
                        patches = varargin{1};
                        n = numel(patches);
                        
                        numbers = getnumber(patches);
                        figure('Name',['Meshes of complementary subdomain and patches #' num2str([numbers{:}])])
                        % set(gcf,'Name',['Meshes of complementary subdomain and patches #' num2str([numbers{:}])])
                        clf
                        hg = cell(1,1+n);
                        leg = cell(1,1+n);
                        % S_out = getmodelpart(glob.S,0);
                        h = plot(glob.S_out,'FaceColor',getfacecolor(1),'LineWidth',p.Results.LineWidth,varargin{:});
                        hg{1} = hggroup;
                        set(h(:),'Parent',hg{1});
                        leg{1} = '$\Omega \setminus \Lambda$';
                        for k=1:n
                            h = plot(patches.patches{k}.S,'FaceColor',getfacecolor(patches.patches{k}.number+1),'LineWidth',p.Results.LineWidth,varargin{:});
                            hg{k+1} = hggroup;
                            set(h(:),'Parent',hg{k+1});
                            if n==1
                                leg{k+1} = '$\Lambda$';
                            else
                                leg{k+1} = ['$\Lambda_{' num2str(patches.patches{k}.number) '}$'];
                            end
                        end
                        if p.Results.legend
                            legend([hg{:}],leg{:},'Location','NorthEastOutside','Interpreter',p.Results.Interpreter)
                        end
                        set(gca,'FontSize',p.Results.FontSize)
                    elseif isa(varargin{1},'Interface')
                        interface = varargin{1};
                        
                        figure('Name',['Meshes of complementary subdomain and interface #' num2str(interface.number)])
                        % set(gcf,'Name',['Meshes of complementary subdomain and interface #' num2str(interface.number)])
                        clf
                        % S_out = getmodelpart(glob.S,0);
                        h_out = plot(glob.S_out,'FaceColor',getfacecolor(1),'LineWidth',p.Results.LineWidth,varargin{:});
                        h_interface = plot(interface.S,'FaceColor',getfacecolor(interface.number+1),'LineWidth',p.Results.LineWidth,varargin{:});
                        hg_out = hggroup;
                        hg_interface = hggroup;
                        set(h_out(:),'Parent',hg_out);
                        set(h_interface(:),'Parent',hg_interface);
                        if p.Results.legend
                            legend([h_out,h_interface],'$\Omega \setminus \Lambda$',['$\Gamma_{' num2str(interface.number) '}$'],...
                                'Location','NorthEastOutside','Interpreter',p.Results.Interpreter)
                        end
                        set(gca,'FontSize',p.Results.FontSize)
                    elseif isa(varargin{1},'Interfaces')
                        interfaces = varargin{1};
                        n = numel(interfaces);
                        
                        numbers = getnumber(interfaces);
                        figure('Name',['Meshes of complementary subdomain and interfaces #' num2str([numbers{:}])])
                        % set(gcf,'Name',['Meshes of complementary subdomain and interfaces #' num2str([numbers{:}])])
                        clf
                        hg = cell(1,1+n);
                        leg = cell(1,1+n);
                        % S_out = getmodelpart(glob.S,0);
                        h = plot(glob.S_out,'FaceColor',getfacecolor(1),'LineWidth',p.Results.LineWidth,varargin{:});
                        hg{1} = hggroup;
                        set(h(:),'Parent',hg{1});
                        leg{1} = '$\Omega \setminus \Lambda$';
                        for k=1:n
                            h = plot(interfaces.interfaces{k}.S,'FaceColor',getfacecolor(interfaces.interfaces{k}.number+1),'LineWidth',p.Results.LineWidth,varargin{:});
                            hg{k+1} = hggroup;
                            set(h(:),'Parent',hg{k+1});
                            if n==1
                                leg{k+1} = '$\Gamma$';
                            else
                                leg{k+1} = ['$\Gamma_{' num2str(interfaces.interfaces{k}.number) '}$'];
                            end
                        end
                        if p.Results.legend
                            legend([hg{:}],leg{:},'Location','NorthEastOutside','Interpreter',p.Results.Interpreter)
                        end
                        set(gca,'FontSize',p.Results.FontSize)
                    end
            end
        end
        
        function plotPartition(glob,varargin)
            % function plotPartition(glob)
            % Display partition of model glob.S
            
            p = ImprovedInputParser;
            addParameter(p,'legend',true,@islogical);
            addParameter(p,'FontSize',16,@isscalar);
            parse(p,varargin{:})
            
            figure('Name','Mesh partition of fictitious domain')
            % set(gcf,'Name','Mesh partition of fictitious domain')
            clf
            plotparamelem(glob.S,'partition');
            if ~p.Results.legend
                legend('off')
            end
            set(gca,'FontSize',p.Results.FontSize)
            
        end
        
    end
    
end
