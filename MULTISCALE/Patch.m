classdef Patch
    %
    
    %
    %   A MATLAB class for representing a patch (local subdomain of
    %   interest)
    %
    %   See also PATCHES, INTERFACE, INTERFACES
    
    properties
        S
        A
        b
        Atang
        M
        C
        b0
        u0
        v0
        number
        type
        changeOfVariable
        increment
        solver
        initializationType
        timeSolver
        timeOrder
        display
    end
    
    methods
        
        % Constructor
        function patch = Patch(varargin)
            % Class Patch
            %
            % function patch = Patch(varargin)
            % patch.S: model
            % patch.A: stiffness operator
            % patch.b: sollicitation vector
            % patch.Atang: tangent stiffness operator
            % patch.M: mass operator
            % patch.C: damping operator
            % patch.b0: initial sollicitation vector, [] by default
            % patch.u0: initial solution, [] by default
            % patch.v0: initial velocity, [] by default
            % patch.number: patch number, [] by default
            % patch.type: patch type, '' by default
            % patch.changeOfVariable: reformulation of local problem as a single-field problem by introducing a change of variable w=w_tilde+z (true or false), false by default
            % patch.increment: reformulation of local problem on increments (instead of current iterates) (true or false), true by default
            % patch.solver: solver for local problem ([], NEWTONSOLVER), [] by default
            % if patch.solver is an instance of class NEWTONSOLVER,
            % patch.initializationType: type of initialization for iterative resolution of local problem ('zero', 'one' or 'previter'), 'zero' by default value
            % patch.timeSolver: time solver for local problem ([], DGTIMESOLVER, EULERSOLVER, NEWMARKSOLVER), [] by default
            % patch.timeOrder: time order for local problem, 0 by default
            % patch.display: display (true or false), false by default
            
            switch nargin
                case 0
                    patch.S = [];
                    patch.A = [];
                    patch.Atang = [];
                    patch.M = [];
                    patch.b = [];
                    patch.b0 = [];
                    patch.u0 = [];
                    patch.v0 = [];
                case 1
                    if isa(varargin{1},'Patch')
                        patch = varargin{1};
                    else
                        patch.S = varargin{1};
                    end
                otherwise
                    patch.S = varargin{1};
                    patch.A = varargin{2};
                    patch.b = varargin{3};
                    patch.Atang = varargin{4};
                    patch.M = varargin{5};
                    patch.C = varargin{6};
                    patch.b0 = varargin{7};
                    patch.u0 = varargin{8};
                    patch.v0 = varargin{9};
            end
            
            expectedInitializationTypes = {'zero','one','previter'};
            
            p = ImprovedInputParser;
            addParameter(p,'number',[],@isscalar);
            addParameter(p,'type','',@ischar);
            addParameter(p,'changeOfVariable',false,@islogical);
            addParameter(p,'increment',true,@islogical);
            addParameter(p,'solver',[]);
            % addParameter(p,'initializationType','zero',@ischar);
            addParameter(p,'initializationType','zero',...
                @(x) any(validatestring(x,expectedInitializationTypes)))
            addParameter(p,'timeSolver',[]);
            addParameter(p,'timeOrder',0);
            addParameter(p,'display',false,@islogical);
            parse(p,varargin{:});
            patch = passMatchedArgsToProperties(p,patch);
            
        end
        
        function number = getnumber(patch)
            % function number = getnumber(patch)
            
            number = patch.number;
        end
        
        function patch = setnumber(patch,number)
            % function patch = setnumber(patch,number)
            
            if ~isa(number,'double')
                error('The number must be a double.')
            end
            patch.number = number;
        end
        
        function [rep,loc] = ismember(k,patch)
            % function [rep,loc] = ismember(k,patch)
            
            if ~isa(k,'double')
                k = getnumber(k);
                if isa(k,'cell')
                    k = [k{:}];
                end
            end
            if ~isa(patch,'double')
                patch = getnumber(patch);
                if isa(patch,'cell')
                    patch = [patch{:}];
                end
            end
            [rep,loc] = ismember(k,patch);
        end
        
        function type = gettype(patch)
            % function type = gettype(patch)
            
            type = patch.type;
        end
        
        function patch = horzcat(varargin)
            % function patch = horzcat(varargin)
            
            patch = Patches(varargin{:});
        end
        
        function [w,lambda,output,vw,vlambda,aw,alambda] = solve(patch,interface,U,w_old,lambda_old,vU,vw_old,vlambda_old,aU,aw_old,alambda_old)
            % function [w,lambda,output] = solve(patch,interface,U,w_old,lambda_old)
            % function [w,lambda,output,vw,vlambda] = solve(patch,interface,U,w_old,lambda_old,vU,vw_old,vlambda_old)
            % function [w,lambda,output,vw,vlambda,aw,alambda] = solve(patch,interface,U,w_old,lambda_old,vU,vw_old,vlambda_old,aU,aw_old,alambda_old)
            % function u = solve(patch,interface,U,w_old,lambda_old)
            % function u = solve(patch,interface,U,w_old,lambda_old,vU,vw_old)
            % function u = solve(patch,interface,U,w_old,lambda_old,vU,vw_old,aU,aw_old)
            % Solves local problem defined over patch
            % If the number of output arguments is equal to one, it returns
            % 1-by-2 cell u containing both current local solution w  and
            % current Lagrange multiplier lambda 
            %
            % Inputs:
            % patch: Patch
            % interface: Interface
            % U: m_U-by-p double containing current global iterate U_{k}
            % w_old: m_w-by-p double containing previous local solution w_{k-1}
            % lambda_old: m_l-by-p double containing previous Lagrange multiplier lambda_{k-1}
            % For first- and second-order time-dependent problems,
            % vU: m_U-by-p double containing current global iterate velocity vU_{k}
            % vw_old: m_w-by-p double containing previous local solution velocity vw_{k-1}
            % vlambda_old: m_l-by-p double containing previous Lagrange multiplier velocity vlambda_{k-1}
            % For second-order time-dependent problems,
            % aU: m_U-by-p double containing current global iterate acceleration aU_{k}
            % aw_old: m_w-by-p double containing previous local solution acceleration aw_{k-1}
            % alambda_old: m_l-by-p double containing previous Lagrange multiplier acceleration alambda_{k-1}
            %
            % Outputs:
            % w: m_w-by-p double containing current local solution w_{k}
            % lambda: m_l-by-p double containing current Lagrange multiplier lambda_{k}
            % u: 1-by-2 cell of double arrays containing current local solution w_{k} and current Lagrange multiplier lambda_{k}
            % output.time : 1-by-1 double containing CPU time t
            % For time-dependent problems,
            % output.result: structure containing outputs of time solver
            % For first- and second-order time-dependent problems,
            % vw: m_w-by-p double containing current local velocity vw_{k}
            % vlambda: m_l-by-p double containing current Lagrange multiplier vlambda_{k}
            % v: 1-by-2 cell containing current local velocity vw_{k} and current Lagrange multiplier vlambda_{k}
            % For second-order time-dependent problems,
            % aw: m_w-by-p double containing current local acceleration aw_{k}
            % alambda: m_l-by-p double containing current Lagrange multiplier alambda_{k}
            % a: 1-by-2 cell containing current local acceleration aw_{k} and current Lagrange multiplier alambda_{k}
            % where
            % k is the current iteration number
            % m_U is the dimension of the spatial approximation space of global solution U
            % m_w is the dimension of the spatial approximation space of local solution w
            % m_l is the dimension of the spatial approximation space of Lagrange multiplier lambda
            % p is the dimension of the time approximation space of global solution U, local solution w and Lagrange multiplier lambda
            
            t = tic;
            
            if nargin<4 || isempty(w_old) || ((nargin<5 || isempty(lambda_old)) && ~patch.changeOfVariable)
                patch.increment = false;
            elseif ~isempty(patch.timeSolver)
                if (nargin<7 || isempty(vw_old) || ((nargin<8 || isempty(vlambda_old)) && ~patch.changeOfVariable)) && patch.timeOrder>=1
                    patch.increment = false;
                end
                if (nargin<10 || isempty(aw_old) || ((nargin<11 || isempty(alambda_old)) && ~patch.changeOfVariable)) && patch.timeOrder>=2
                    patch.increment = false;
                end
            end
            
            patch = initializeRightHandSide(patch);
            
            if patch.changeOfVariable
                patch = setBoundaryCondition(patch,interface);
            end
            
            % Stiffness matrix A
            [A,Atang] = assembleStiffnessOperator(patch,interface);
            
            if ~isempty(patch.timeSolver)
                % Mass matrix M
                M = assembleMassOperator(patch,interface);
                if patch.timeOrder==2
                    % Damping matrix C
                    C = assembleDampingOperator(patch,interface);
                end
            end
            
            % Sollicitation vector b
            b = assembleRightHandSide(patch,interface,U);
            
            if patch.increment && (isempty(patch.solver) || ~isanlsolver(patch.solver))
                u_old = setSolution(patch,w_old,lambda_old,interface,U);
                b = b - A*u_old;
                if ~isempty(patch.timeSolver)
                    v_old = setSolution(patch,vw_old,vlambda_old,interface,vU);
                    if patch.timeOrder==1
                        b = b - M*v_old;
                    elseif patch.timeOrder==2
                        a_old = setSolution(patch,aw_old,alambda_old,interface,aU);
                        b = b - C*v_old - M*a_old;
                    end
                end
            end
            
            % Solution u
            if isempty(patch.timeSolver)
                % Time-independent problem
                if isempty(patch.solver)
                    u = A\b;
                else
                    if patch.changeOfVariable && isanlsolver(patch.solver)
                        warning('Wrong local transfert from patch to interface. Set parameter ''changeOfVariable'' to false.')
                    end
                    if isa(patch.solver,'NEWTONSOLVER')
                        u_old = setSolution(patch,w_old,lambda_old,interface,U);
                        u0 = setInitialSolution(patch.initializationType,b,u_old);
                        u = solve(patch.solver,b,A,Atang,u0);
                    end
                end
            else
                % Time-dependent problem
                if patch.timeOrder==1
                    u0 = setInitialCondition(patch,interface);
                    [u,output.result,v] = dsolve(patch.timeSolver,b,M,A,u0);
                elseif patch.timeOrder==2
                    [u0,v0] = setInitialCondition(patch,interface);
                    [u,output.result,v,a] = ddsolve(patch.timeSolver,b,M,A,C,u0,v0);
                end
            end
            
            if patch.increment && (isempty(patch.solver) || ~isanlsolver(patch.solver))
                u = u_old + u;
                if ~isempty(patch.timeSolver)
                    if patch.timeOrder>=1
                        v = v_old + v;
                    end
                    if patch.timeOrder>=2
                        a = a_old + a;
                    end
                end
            end
            
            if ~patch.changeOfVariable
                % Without change of variable
                % Local solution (w,lambda)
                [w,lambda] = getSolution(patch,u);
                if ~isempty(patch.timeSolver)
                    if patch.timeOrder>=1
                        [vw,vlambda] = getSolution(patch,v);
                    end
                    if patch.timeOrder>=2
                        [aw,alambda] = getSolution(patch,a);
                    end
                end
            else
                % With change of variable
                % Local solution w
                % Lagrange multiplier lambda
                w = getSolution(patch,u,interface,U);
                if isempty(patch.timeSolver)
                    lambda = computeLagrangeMultiplier(patch,interface,w);
                else
                    vw = getSolution(patch,v,interface,vU);
                    if patch.timeOrder==1
                        lambda = computeLagrangeMultiplier(patch,interface,w,vw);
                    elseif patch.timeOrder==2
                        aw = getSolution(patch,a,interface,aU);
                        lambda = computeLagrangeMultiplier(patch,interface,w,vw,aw);
                    end
                    T = gettimemodel(patch.timeSolver);
                    tSolver = EULERTIMESOLVER(T,'eulertype','implicit','display',false);
                    if patch.timeOrder>=1
                        vlambda = derivative(tSolver,lambda);
                    end
                    if patch.timeOrder>=2
                        alambda = derivative(tSolver,vlambda);
                    end
                end
                
                u = vertcat(w,lambda);
                if ~isempty(patch.timeSolver)
                    if patch.timeOrder>=1
                        v = vertcat(vw,vlambda);
                    end
                    if patch.timeOrder>=2
                        a = vertcat(aw,alambda);
                    end
                end
            end
            
            output.time = toc(t);
            
            if patch.display
                fprintf('\nElapsed time = %f s\n',output.time);
            end
            
            if nargout==1
                w = cellOutput(patch,interface,u);
            end
            
        end
        
        function patch = getPatch(patch,k)
            % function patch = getPatch(patch,k)
            
            if nargin==2 && k~=patch.number
                error('Wrong patch number.')
            end
        end
        
        function patch = eval(patch,xi,varargin)
            % function patch = eval(patch,xi,varargin)
            
            if israndom(patch.S)
                patch.S = randomeval(patch.S,xi,varargin{:});
            end
            if israndom(patch.A)
                patch.A = randomeval(patch.A,xi,varargin{:});
            end
            if israndom(patch.b)
                patch.b = randomeval(patch.b,xi,varargin{:});
            end
            if israndom(patch.Atang)
                patch.Atang = randomeval(patch.Atang,xi,varargin{:});
            end
            if israndom(patch.M)
                patch.M = randomeval(patch.M,xi,varargin{:});
            end
            if israndom(patch.C)
                patch.C = randomeval(patch.C,xi,varargin{:});
            end
            if israndom(patch.b0)
                patch.b0 = randomeval(patch.b0,xi,varargin{:});
            end
            if israndom(patch.u0)
                patch.u0 = randomeval(patch.u0,xi,varargin{:});
            end
            if israndom(patch.v0)
                patch.v0 = randomeval(patch.v0,xi,varargin{:});
            end
        end

        function patch = calcOperator(patch,varargin)
            % function patch = calcOperator(patch,varargin)
            
            if ~isempty(patch.timeSolver)
                patch.M = calc_mass(patch.S);
                if ~isempty(patch.C)
                    patch.C = calc_damp(patch.C);
                end
            end
            if isempty(patch.solver)
                if ~length(MATERIALS(patch.S))
                    a_patch = varargin{1};
                    patch.A = calc_matrix(a_patch{patch.number},patch.S);
                else
                    patch.A = calc_rigi(patch.S);
                end
            else
                if isa(patch.solver,'NEWTONSOLVER')
                    patch.A = @(u) calc_fint(patch.S,u);
                    patch.Atang = @(u) calc_rigitang(patch.S,u);
                end
            end
        end
        
        function patch = linearizeOperator(patch,u)
            % function patch = linearizeOperator(patch,u)
            
            if ~isempty(patch.solver) && isanlsolver(patch.solver)
                patch.solver = [];
                patch.A = patch.Atang(u);
            end
        end
        
        function patch = initializeRightHandSide(patch)
            % function patch = initializeRightHandSide(patch)
            
            if isempty(patch.b)
                sz = getnbddl(patch.S);
                patch.b = sparse(sz,1);
                if ~isempty(patch.timeSolver)
                    T = gettimemodel(patch.timeSolver);
                    patch.b = patch.b*zero(T);
                end
            end
        end
        
        function b = assembleRightHandSide(patch,interface,U)
            % function b = assembleRightHandSide(patch,interface,U)
            
            if ~patch.changeOfVariable
                % Without change of variable
                B_glob = interface.P_glob'*interface.M;
                
                b1 = patch.b;
                b2 = -B_glob'*U;
                b = vertcat(b1,b2);
            else
                % With change of variable
                P = interface.P_patch'*interface.P_glob;
                
                if isempty(patch.solver)
                    b = freevector(patch.S,patch.b - patch.A*P*U);
                else
                    if isa(patch.solver,'NEWTONSOLVER')
                        b = freevector(patch.S,patch.b - patch.A(P*U));
                    end
                end
                
                if ~isempty(patch.timeSolver)
                    if patch.timeOrder==1
                        b = b - freevector(patch.S,patch.M*P*vU);
                    elseif patch.timeOrder==2
                        b = b - freevector(patch.S,patch.M*P*aU);
                        if ~isempty(patch.C)
                            b = b - freevector(patch.S,patch.C*P*vU);
                        end
                    end
                end
            end
        end
        
        function [A,Atang] = assembleStiffnessOperator(patch,interface)
            % function [A,Atang] = assembleStiffnessOperator(patch,interface)
            
            if ~patch.changeOfVariable
                % Without change of variable
                sz_w = getnbddl(patch.S);
                sz_lambda = getnbddl(interface.S);
                
                B_patch = interface.P_patch'*interface.M;
                
                if isempty(patch.solver)
                    A11 = patch.A;
                    A22 = sparse(sz_lambda,sz_lambda);
                    A12 = -B_patch;
                    A21 = A12';
                    A = [A11,A12; A21,A22];
                    Atang = [];
                else
                    if isa(patch.solver,'NEWTONSOLVER')
                        A1 = @(w,lambda) patch.A(w) - B_patch*lambda;
                        A2 = @(w) -B_patch'*w;
                        A = @(u) [A1(u(1:sz_w),u(sz_w+1:end)); A2(u(1:sz_w))];
                        
                        Atang11 = @(w) patch.Atang(w);
                        Atang22 = sparse(sz_lambda,sz_lambda);
                        Atang12 = -B_patch;
                        Atang21 = Atang12';
                        Atang = @(u) [Atang11(u(1:sz_w)),Atang12; Atang21,Atang22];
                    end
                end
            else
                % With change of variable
                if isempty(patch.solver)
                    A = freematrix(patch.S,patch.A);
                    Atang = [];
                else
                    if isa(patch.solver,'NEWTONSOLVER')
                        A = @(z) freevector(patch.S,patch.A(z));
                        Atang = @(z) freematrix(patch.S,patch.Atang(z));
                    end
                end
            end
        end
        
        function M = assembleMassOperator(patch,interface)
            % function M = assembleMassOperator(patch,interface)
            
            if ~patch.changeOfVariable
                % Without change of variable
                sz_w = getnbddl(patch.S);
                sz_lambda = getnbddl(interface.S);
                
                M11 = patch.M;
                M22 = sparse(sz_lambda,sz_lambda);
                M12 = sparse(sz_w,sz_lambda);
                M21 = M12';
                M = [M11,M12; M21,M22];
            else
                % With change of variable
                M = freematrix(patch.S,patch.M);
            end
        end
        
        function C = assembleDampingOperator(patch,interface)
            % function C = assembleDampingOperator(patch,interface)
            
            if ~patch.changeOfVariable
                % Without change of variable
                sz_w = getnbddl(patch.S);
                sz_lambda = getnbddl(interface.S);
                
                if ~isempty(patch.C)
                    C11 = patch.C;
                else
                    C11 = sparse(sz_w,sz_w);
                end
                C22 = sparse(sz_lambda,sz_lambda);
                C12 = sparse(sz_w,sz_lambda);
                C21 = C12';
                C = [C11,C12; C21,C22];
            else
                % With change of variable
                if ~isempty(patch.C)
                    C = freematrix(patch.S,patch.C);
                else
                    sz = getnbddlfree(patch.S);
                    C = sparse(sz,sz);
                end
            end
        end
        
        function patch = setBoundaryCondition(patch,interface)
            % function patch = setBoundaryCondition(patch,interface)
            
            [~,numnode,~] = intersect(patch.S,interface.S);
            patch.S = addcl(patch.S,numnode);  
        end
        
        function [u0,v0] = setInitialCondition(patch,interface)
            % function u0 = setInitialCondition(patch,interface)
            % function [u0,v0] = setInitialCondition(patch,interface)
            
            if ~patch.changeOfVariable
                % Without change of variable
                % Initial local solution u0=(w0,lambda0)
                w0 = initializeSolution(patch);
                lambda0 = initializeSolution(interface);
                u0 = [w0;lambda0];
                if nargout>1
                    % Initial local velocity v0=(vw0,vlambda0)
                    vw0 = initializeVelocity(patch);
                    vlambda0 = initializeVelocity(interface);
                    v0 = [vw0;vlambda0];
                end
            else
                % With change of variable
                sz = getnbddlfree(patch.S);
                % Initial local solution u0
                if isempty(patch.u0)
                    u0 = zeros(sz,1);
                else
                    u0 = freevector(patch.S,patch.u0);
                end
                if nargout>1
                    % Initial local velocity v0
                    if isempty(patch.v0)
                        v0 = zeros(sz,1);
                    else
                        v0 = freevector(patch.S,patch.v0);
                    end
                end
            end
        end
        
        function w0 = initializeSolution(patch)
            % function w0 = initializeSolution(patch)
            
            sz_w = getnbddl(patch.S);
            if isempty(patch.u0)
                w0 = zeros(sz_w,1);
            else
                w0 = patch.u0;
            end
        end
        
        function vw0 = initializeVelocity(patch)
            % function vw0 = initializeVelocity(interface)
            
            sz = getnbddl(patch.S);
            if isempty(patch.v0)
                vw0 = zeros(sz,1);
            else
                vw0 = patch.v0;
            end
        end
        
        function u = setSolution(patch,w,lambda,interface,U)
            % function u = setSolution(patch,w,lambda)
            % function u = setSolution(patch,w,lambda,interface,U)
            
            if ~patch.changeOfVariable
                % Without change of variable
                u = vertcat(w,lambda);
            else
                % With change of variable
                P = interface.P_patch'*interface.P_glob;
                u = freevector(patch.S,w - P*U);
            end
        end
        
        function [w,lambda] = getSolution(patch,u,interface,U)
            % function [w,lambda] = getSolution(patch,u,interface,U)
            % function w = getSolution(patch,u,interface,U)
            
            if ~patch.changeOfVariable
                % Without change of variable
                sz = getnbddl(patch.S);
                w = u(1:sz);
                lambda = u(sz+1:end);
            else
                % With change of variable
                P = interface.P_patch'*interface.P_glob;
                w = unfreevector(patch.S,u) + P*U;
            end
        end
        
        function lambda = computeLagrangeMultiplier(patch,interface,w,vw,aw)
            % function lambda = computeLagrangeMultiplier(patch,interface,w)
            % function lambda = computeLagrangeMultiplier(patch,interface,w,vw)
            % function lambda = computeLagrangeMultiplier(patch,interface,w,vw,aw)
            
            % Stiffness matrix A
            % B_patch = interface.P_patch'*interface.M;
            % A = B_patch;
            A = interface.M;
            
            % Sollicitation vector b
            if isempty(patch.solver) || ~isanlsolver(patch.solver)
                b = patch.A*w - patch.b;
            else
                b = patch.A(w) - patch.b;
            end
            if ~isempty(patch.timeSolver)
                if patch.timeOrder==1
                    b = b + patch.M*vw;
                elseif patch.timeOrder==2
                    b = b + patch.M*aw;
                    if ~isempty(patch.C)
                        b = b + patch.C*vw;
                    end
                end
            end
            b = interface.P_patch*b;
            
            % Lagrange multiplier lambda
            if isempty(patch.timeSolver)
                lambda = A\b;
            else
                T = gettimemodel(w);
                b_val = getvalue(b);
                lambda = A\b_val;
                lambda = TIMEMATRIX(lambda,T);
            end
        end
        
        function u = cellOutput(patch,interface,u)
            % function u = cellOutput(patch,interface,u)
            
            sz_w = getnbddl(patch.S);
            sz_lambda = getnbddl(interface.S);
            rep = [sz_w,sz_lambda];
            if isempty(patch.timeSolver)
                u = mat2cell(u,rep);
            else
                % T = gettimemodel(u);
                u_val = getvalue(u);
                u = mat2cell(u_val,rep);
                % u = cell(1,2);
                % u_val = mat2cell(u_val,rep);
                % u{1} = TIMEMATRIX(u_val{1},T);
                % u{2} = TIMEMATRIX(u_val{2},T);
            end
        end
        
        function [P,numnode] = calcProjection(patch,interface)
            % function [P,numnode] = calcProjection(patch,interface)
            % Calculates projection operator from model patch.S to model interface.S
            % patch: Patch
            % interface: Interface
            
            [P,numnode] = calc_P_free(patch.S,interface.S);
        end
        
        function plotModel(patch,varargin)
            % function plotModel(patch)
            % Display model patch.S
            
            p = ImprovedInputParser;
            addParameter(p,'legend',true,@islogical);
            addParameter(p,'FontSize',16,@isscalar);
            addParameter(p,'LineWidth',0.5,@isscalar);
            addParameter(p,'Interpreter','latex',@ischar);
            parse(p,varargin{:})
            
            varargin = delcharin({'legend','FontSize','LineWidth','Interpreter'},varargin);
            
            figure('Name',['Mesh of patch #' num2str(patch.number)])
            % set(gcf,'Name',['Mesh of patch #' num2str(patch.number)])
            clf
            h = plot(patch.S,'FaceColor',getfacecolor(patch.number+1),'LineWidth',p.Results.LineWidth,varargin{:});
            hg = hggroup;
            set(h(:),'Parent',hg);
            if p.Results.legend
                legend(hg,['$\Lambda_{' num2str(patch.number) '}$'],'Location','NorthEastOutside','Interpreter',p.Results.Interpreter)
            end
            set(gca,'FontSize',p.Results.FontSize)
            
        end

    end
    
    methods (Static)
        
        function n = numel(~)
            % function n = numel(patch)
            n = 1;
        end
        
    end
    
end
