classdef GlobalOutside
    %
    
    %
    %   A MATLAB class for representing the global problem defined over the
    %   complementary subdomain (global domain outside patches)
    %
    
    properties( SetAccess = public, GetAccess = public )
        S
        A
        b
        Atang
        M
        C
        b0
        u0
        v0
        increment
        solver
        initializationType
        timeSolver
        timeOrder
        display
    end
    
    methods( Access = public )
        
        % Constructor
        function glob = GlobalOutside(varargin)
            % Class GlobalOutside
            %
            % function glob = GlobalOutside(varargin)
            % glob.S: model associated to complementary subdomain
            % glob.A: stiffness operator associated to complementary subdomain
            % glob.b: sollicitation vector associated to complementary subdomain
            % glob.Atang: tangent stiffness operator associated to complementary subdomain
            % glob.M: mass operator associated to complementary subdomain
            % glob.C: damping operator associated to complementary subdomain
            % glob.b0: initial sollicitation vector associated to complementary subdomain, [] by default
            % glob.u0: initial solution associated to complementary subdomain, [] by default
            % glob.v0: initial velocity associated to complementary subdomain, [] by default
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
                    glob.A = [];
                    glob.b = [];
                    glob.Atang = [];
                    glob.M = [];
                    glob.C = [];
                    glob.b0 = [];
                    glob.u0 = [];
                    glob.v0 = [];
                case 1
                    if isa(varargin{1},'GlobalOutside')
                        glob = varargin{1};
                    else
                        glob.S = varargin{1};
                    end
                otherwise
                    glob.S = varargin{1};
                    glob.A = varargin{2};
                    glob.b = varargin{3};
                    glob.Atang = varargin{4};
                    glob.M = varargin{5};
                    glob.C = varargin{6};
                    glob.b0 = varargin{7};
                    glob.u0 = varargin{8};
                    glob.v0 = varargin{9};
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
        
        function [V,output,vV,aV] = solve(glob,interfaces,lambda_old,V_old,vV_old,aV_old)
            % function [V,output] = solve(glob,interfaces,lambda_old,V_old)
            % function [V,output,vV] = solve(glob,interfaces,lambda_old,V_old,vV_old)
            % function [V,output,vV,aV] = solve(glob,interfaces,lambda_old,V_old,vV_old,aV_old)
            % Solves global problem defined over complementary subdomain
            %
            % Inputs:
            % glob: GlobalOutside
            % interfaces: Interfaces
            % lambda_old: 1-by-n cell of m_l-by-p double containing previous Lagrange multiplier lambda_{k-1}
            % V_old: m_U-by-p double containing previous global solution V_{k-1}
            % For first- and second-order time-dependent problems,
            % vV_old: m_U-by-p double containing previous global solution velocity vV_{k-1}
            % For second-order time-dependent problems,
            % aV_old: m_U-by-p double containing previous global solution acceleration aV_{k-1}
            %
            % Outputs:
            % V: m_U-by-p double containing current global solution V_{k}
            % output.time : 1-by-1 double containing CPU time t
            % For time-dependent problems,
            % output.result: structure containing outputs of time solver
            % For first- and second-order time-dependent problems,
            % vV: m_U-by-p double containing current global velocity vV_{k}
            % For second-order time-dependent problems,
            % aV: m_U-by-p double containing current global acceleration aV_{k}
            % where
            % k is the current iteration number
            % m_U is the dimension of the spatial approximation space of global solution U
            % m_l is the dimension of the spatial approximation space of Lagrange multiplier lambda
            % p is the dimension of the time approximation space of global solution U and Lagrange multiplier lambda
            
            t = tic;
            
            if nargin<4 || isempty(V_old)
                glob.increment = false;
            elseif ~isempty(glob.timeSolver)
                if (nargin<6 || isempty(vV_old)) && glob.timeOrder>=1
                    glob.increment = false;
                end
                if (nargin<8 || isempty(aV_old)) && glob.timeOrder>=2
                    glob.increment = false;
                end
            end
            
            glob = initializeRightHandSide(glob);
            b = setRightHandSide(glob,interfaces,lambda_old);
            
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
            % glob: GlobalOutside
            % U: global solution
            
            if ~isempty(glob.solver) && isanlsolver(glob.solver)
                glob.solver = [];
                glob.A = glob.Atang(U);
            end
            
        end
        
        function glob = initializeRightHandSide(glob)
            % function glob = initializeRightHandSide(glob)
            % Initialize right-hand side of global problem
            
            if isempty(glob.b)
                sz = getnbddlfree(glob.S);
                glob.b = sparse(sz,1);
                if ~isempty(glob.timeSolver)
                    T = gettimemodel(glob.timeSolver);
                    glob.b = glob.b*zero(T);
                end
            end
        end
        
        function b = setRightHandSide(glob,interfaces,lambda_old)
            % function b = setRightHandSide(glob,interfaces,lambda_old)
            % Set right-hand side of global problem
            % glob: GlobalOutside
            % interfaces: Interfaces
            % lambda_old: 1-by-n cell of m_l-by-p double containing previous Lagrange multiplier lambda_{k-1}
            
            b = glob.b;
            n = numel(interfaces);
            interface = interfaces.interfaces;
            for k=1:n
                B_glob = interface{k}.P_globOut'*interface{k}.M;
                b = b - B_glob*lambda_old{k};
            end
        end
        
        function [U0,V0] = setInitialCondition(glob)
            % function U0 = setInitialCondition(glob)
            % Sets global initial solution vector U0
            % function [U0,V0] = setInitialCondition(glob)
            % Sets global initial solution vector U0 and velocity vector V0
            % glob: GlobalOutside
            
            sz = getnbddlfree(glob.S);
            if isempty(glob.u0)
                U0 = zeros(sz,1);
            else
                U0 = glob.u0;
            end
            if nargout>1
                if isempty(glob.v0)
                    V0 = zeros(sz,1);
                else
                    V0 = glob.v0;
                end
            end
        end
        
        function [P,numnode] = calcProjection(glob,varargin)
            % function [P,numnode] = calcProjection(glob,varargin)
            % Calculates projection operator from model glob.S to model varargin{1}.S if varargin{1} is an instance of type Interface
            % Calculates projection operator from model glob.S to model varargin{1} if varargin{1} is an instance of class MODEL
            % glob: GlobalOutside
            
            if isa(varargin{1},'Interface')
                S_interface = varargin{1}.S;
                [P,numnode] = calcProjection(S_interface,glob.S,varargin{2:end});
                P = P';
            elseif isa(varargin{1},'MODEL')
                free = getcharin('free',varargin,true);
                varargin = delcharin('free',varargin);
                if free
                    [P,numnode] = calc_P_free(glob.S,varargin{:});
                else
                    [P,numnode] = calc_P(glob.S,varargin{:});
                end
            end
            
        end
        
        function plotModel(glob,varargin)
            % function plotModel(glob)
            % Display model glob.S
            %
            % function plotModel(glob,patch)
            % Display model glob.S as well as model patch.S
            %
            % function plotModel(glob,patches)
            % Display model glob.S as well as models patch.S of all patches
            %
            % function plotModel(glob,interface)
            % Display model glob.S as well as model interface.S
            %
            % function plotModel(glob,interfaces)
            % Display model glob.S as well as models interface.S of all interfaces
            
            p = ImprovedInputParser;
            addParameter(p,'legend',true,@islogical);
            addParameter(p,'FontSize',16,@isscalar);
            addParameter(p,'LineWidth',0.5,@isscalar);
            addParameter(p,'Interpreter','latex',@ischar);
            if isempty(varargin) || ischar(varargin{1})
                parse(p,varargin{:})
            else
                parse(p,varargin{2:end})
            end
            
            varargin = delcharin({'legend','FontSize','LineWidth','Interpreter'},varargin);
            nargin = length(varargin)+1;
            
            switch nargin
                case 1
                    figure('Name','Mesh of complementary subdomain')
                    % set(gcf,'Name','Mesh of complementary subdomain')
                    clf
                    h = plot(glob.S,'FaceColor',getfacecolor(1),'LineWidth',p.Results.LineWidth,varargin{:});
                    hg = hggroup;
                    set(h(:),'Parent',hg);
                    if p.Results.legend
                        legend(hg,'$\Omega \setminus \Lambda$','Location','NorthEastOutside','Interpreter',p.Results.Interpreter)
                    end
                    set(gca,'FontSize',p.Results.FontSize)
                case 2
                    if isa(varargin{1},'Patch')
                        patch = varargin{1};
                        
                        figure('Name',['Meshes of complementary subdomain and patch #' num2str(patch.number)])
                        % set(gcf,'Name',['Meshes of complementary subdomain and patch #' num2str(patch.number)])
                        clf
                        h = plot(glob.S,'FaceColor',getfacecolor(1),'LineWidth',p.Results.LineWidth,varargin{:});
                        h_patch = plot(patch.S,'FaceColor',getfacecolor(patch.number+1),'LineWidth',p.Results.LineWidth,varargin{:});
                        hg = hggroup;
                        hg_patch = hggroup;
                        set(h(:),'Parent',hg);
                        set(h_patch(:),'Parent',hg_patch);
                        if p.Results.legend
                            legend([hg,hg_patch],'$\Omega \setminus \Lambda$',['$\Lambda_{' num2str(patch.number) '}$'],...
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
                        h = plot(glob.S,'FaceColor',getfacecolor(1),'LineWidth',p.Results.LineWidth,varargin{:});
                        hg{1} = hggroup;
                        set(h(:),'Parent',hg{1});
                        leg{1} = '$\Omega \setminus \Lambda$';
                        for k=1:n
                            patch = patches.patches{k};
                            h = plot(patch.S,'FaceColor',getfacecolor(patch.number+1),'LineWidth',p.Results.LineWidth,varargin{:});
                            hg{k+1} = hggroup;
                            set(h(:),'Parent',hg{k+1});
                            if n==1
                                leg{k+1} = '$\Lambda$';
                            else
                                leg{k+1} = ['$\Lambda_{' num2str(patch.number) '}$'];
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
                        h = plot(glob.S,'FaceColor',getfacecolor(1),'LineWidth',p.Results.LineWidth,varargin{:});
                        h_interface = plot(interface.S,'FaceColor',getfacecolor(interface.number+1),'LineWidth',p.Results.LineWidth,varargin{:});
                        hg = hggroup;
                        hg_interface = hggroup;
                        set(h(:),'Parent',hg);
                        set(h_interface(:),'Parent',hg_interface);
                        if p.Results.legend
                            legend([hg,hg_interface],'$\Omega \setminus \Lambda$',['$\Gamma_{' num2str(interface.number) '}$'],...
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
                        h = plot(glob.S,'FaceColor',getfacecolor(1),'LineWidth',p.Results.LineWidth,varargin{:});
                        hg{1} = hggroup;
                        set(h(:),'Parent',hg{1});
                        leg{1} = '$\Omega \setminus \Lambda$';
                        for k=1:n
                            interface = interfaces.interfaces{k};
                            h = plot(interface.S,'FaceColor',getfacecolor(interface.number+1),'LineWidth',p.Results.LineWidth,varargin{:});
                            hg{k+1} = hggroup;
                            set(h(:),'Parent',hg{k+1});
                            if n==1
                                leg{k+1} = '$\Gamma$';
                            else
                                leg{k+1} = ['$\Gamma_{' num2str(interface.number) '}$'];
                            end
                        end
                        if p.Results.legend
                            legend([hg{:}],leg{:},'Location','NorthEastOutside','Interpreter',p.Results.Interpreter)
                        end
                        set(gca,'FontSize',p.Results.FontSize)
                    end
                    
            end
            
        end
        
    end
    
end
