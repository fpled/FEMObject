classdef Interface
    %
    
    %
    %   A MATLAB class for representing an interface
    %
    %   See also INTERFACES, PATCH, PATCHES
    
    properties
        S
        M
        u0
        v0
        P_glob
        P_globOut
        P_patch
        number
        type
    end
    
    methods
        
        % Constructor
        function interface = Interface(varargin)
            % Class Interface
            %
            % function interface = Interface(varargin)
            % interface.S: model
            % interface.M: mass matrix
            % interface.u0: initial solution, [] by default
            % interface.v0: initial velocity, [] by default
            % interface.P_glob: projection operator from ficitious domain to interface
            % interface.P_globOut: projection operator from complementary subdomain to interface
            % interface.P_patch: projection operator from patch to interface
            % interface.number: interface number, [] by default
            % interface.type: interface type, '' by default
            %
            % Example
            % patch = Patch();
            % interface = Interface(patch);
            
            p = ImprovedInputParser;
            addParameter(p,'number',[],@isscalar);
            addParameter(p,'type','',@ischar);
            
            switch nargin
                case 0
                    interface.S = [];
                    interface.M = [];
                    interface.u0 = [];
                    interface.v0 = [];
                    interface.P_glob = [];
                    interface.P_globOut = [];
                    interface.P_patch = [];
                    parse(p,varargin{:});
                case 1
                    if isa(varargin{1},'Interface')
                        interface = varargin{1};
                    elseif isa(varargin{1},'Patch')
                        patch = varargin{1};
                        interface.S = create_boundary(patch.S);
                        interface.S = final(interface.S);
                    else
                        interface.S = varargin{1};
                    end
                    parse(p,varargin{2:end});
                case 2
                    if isa(varargin{1},'Interface')
                        interface = varargin{1};
                    else
                        if isa(varargin{1},'Patch') && (isa(varargin{2},'Global') || isa(varargin{2},'GlobalOutside'))
                            patch = varargin{1};
                            glob = varargin{2};
                        elseif (isa(varargin{1},'Global') || isa(varargin{1},'GlobalOutside')) && isa(varargin{2},'Patch')
                            glob = varargin{1};
                            patch = varargin{2};
                        end
                        interface.S = create_boundary(patch.S);
                        if isa(glob,'Global')
                            interface.S = intersect(interface.S,glob.S_out);
                        elseif isa(glob,'GlobalOutside')
                            interface.S = intersect(interface.S,glob.S);
                        end
                        interface.S = final(interface.S);
                    end
                    parse(p,varargin{3:end});
                otherwise
                    interface.S = varargin{1};
                    interface.M = varargin{2};
                    interface.u0 = varargin{3};
                    interface.v0 = varargin{4};
                    interface.P_glob = varargin{5};
                    interface.P_globOut = varargin{6};
                    interface.P_patch = varargin{7};
                    parse(p,varargin{8:end});
            end

            interface = passMatchedArgsToProperties(p,interface);
            
        end
        
        function number = getnumber(interface)
            % function number = getnumber(interface)
            
            number = interface.number;
        end
        
        function interface = setnumber(interface,number)
            % function interface = setnumber(interface,number)
            
            if ~isa(number,'double')
                error('The number must be a double.')
            end
            interface.number = number;
        end
        
        function [rep,loc] = ismember(k,interface)
            % function [rep,loc] = ismember(k,interface)
            
            if ~isa(k,'double')
                k = getnumber(k);
                if isa(k,'cell')
                    k = [k{:}];
                end
            end
            if ~isa(interface,'double')
                interface = getnumber(interface);
                if isa(interface,'cell')
                    interface = [interface{:}];
                end
            end
            [rep,loc] = ismember(k,interface);
        end
        
        function type = gettype(interface)
            % function type = gettype(interface)
            
            type = interface.type;
        end
        
        function interface = horzcat(varargin)
            % function interface = horzcat(varargin)
            
            interface = Interfaces(varargin{:});
        end
        
        function lambda0 = initializeSolution(interface)
            % function lambda0 = initializeSolution(interface)
            
            sz = getnbddl(interface.S);
            if isempty(interface.u0)
                lambda0 = zeros(sz,1);
            else
                lambda0 = patch.u0;
            end
        end
        
        function vlambda0 = initializeVelocity(interface)
            % function vlambda0 = initializeVelocity(interface)
            
            sz = getnbddl(interface.S);
            if isempty(interface.v0)
                vlambda0 = zeros(sz,1);
            else
                vlambda0 = interface.v0;
            end
        end
        
        function [P,numnode] = calcProjection(interface,glob,varargin)
            % function [P,numnode] = calcProjection(interface,glob,varargin)
            % Calculates projection operator from model interface.S to model glob.S
            % interface: Interface
            % glob: Global or GlobalOutside
            
            [P,numnode] = calcProjection(interface.S,glob.S,varargin{:});
            
        end
        
        function plotModel(interface,varargin)
            % function plotModel(interface)
            % Display model interface.S
            
            p = ImprovedInputParser;
            addParameter(p,'legend',true,@islogical);
            addParameter(p,'FontSize',16,@isscalar);
            addParameter(p,'LineWidth',0.5,@isscalar);
            addParameter(p,'Interpreter','latex',@ischar);
            parse(p,varargin{:})
            
            varargin = delcharin({'legend','FontSize','LineWidth','Interpreter'},varargin);
            
            figure('Name',['Mesh of interface #' num2str(interface.number)])
            % set(gcf,'Name',['Mesh of interface #' num2str(Interface.number)])
            clf
            h = plot(interface.S,'FaceColor',getfacecolor(interface.number+1),'LineWidth',p.Results.LineWidth,varargin{:});
            hg = hggroup;
            set(h(:),'Parent',hg);
            if p.Results.legend
                legend(hg,['$\Gamma_{' num2str(interface.number) '}$'],'Location','NorthEastOutside','Interpreter',p.Results.Interpreter)
            end
            set(gca,'FontSize',p.Results.FontSize)
            
        end

    end
    
    methods (Static)
        
        function n = numel(~)
            % function n = numel(interface)
            n = 1;
        end
        
    end
    
end
