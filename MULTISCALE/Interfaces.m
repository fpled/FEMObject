classdef Interfaces
    %
    
    %
    %   A MATLAB class for representing a collection of interfaces
    %
    %   See also INTERFACE, PATCH, PATCHES
    
    properties
        interfaces
    end
    
    methods
        
        % Constructor
        function interfaces = Interfaces(varargin)
            % Class Interfaces
            %
            % function interfaces = Interfaces(interfaces)
            % interfaces: cell containing objects of type Interface
            %
            % function interfaces = Interfaces(n)
            % interfaces: cell of length n containing objects of type Interface
            
            switch nargin
                case 0
                    interfaces.interfaces = {};
                case 1
                    if isa(varargin{1},'Interfaces')
                        interfaces = varargin{1};
                    elseif isa(varargin{1},'Patches')
                        patches = varargin{1};
                        n = numel(patches);
                        for k=1:n
                            patch = getPatch(patches,k);
                            interfaces.interfaces{k} = Interface(patch);
                        end
                        interfaces = setnumber(interfaces);
                    elseif isa(varargin{1},'double')
                        n = varargin{1};
                        interfaces.interfaces = cell(1,n);
                        for k=1:n
                            interfaces.interfaces{k} = Interface();
                        end
                        interfaces = setnumber(interfaces);
                    end
                case 2
                    if isa(varargin{1},'Interfaces')
                        interfaces = varargin{1};
                    elseif isa(varargin{1},'Patches') && isa(varargin{2},'Global')
                        patches = varargin{1};
                        n = numel(patches);
                        for k=1:n
                            patch = getPatch(patches,k);
                            interfaces.interfaces{k} = Interface(patch,varargin{2});
                        end
                        interfaces = setnumber(interfaces);
                    elseif isa(varargin{1},'double')
                        n = varargin{1};
                        interfaces.interfaces = cell(1,n);
                        for k=1:n
                            interfaces.interfaces{k} = Interface();
                        end
                        interfaces = setnumber(interfaces);
                    end
                otherwise
                    interfaces.interfaces = cell(1,0);
                    for k=1:nargin
                        if isa(varargin{k},'Interface')
                            interfaces.interfaces = [interfaces.interfaces, {varargin{k}}];
                        elseif isa(varargin{k},'Interfaces')
                            interfaces.interfaces = [interfaces.interfaces, varargin{k}.interfaces];
                        elseif isa(varargin{k},'cell')
                            interfaces_k = Interfaces(varargin{k}{:}) ;
                            interfaces.interfaces = [interfaces.interfaces, interfaces_k.interfaces];
                        elseif isa(varargin{k},'Patch')
                            interfaces.interfaces = [interfaces.interfaces, Interface(varargin{k})];
                        elseif isa(varargin{k},'Patches')
                            interfaces.interfaces = [interfaces.interfaces, Interfaces(varargin{k})];
                        end
                    end
                    interfaces = setnumber(interfaces);
                    interfaces = unique(interfaces);
            end
            
        end
        
        function n = numel(interfaces)
            % function n = numel(interfaces)
            
            n = numel(interfaces.interfaces);
        end
        
        function number = getnumber(interfaces,k)
            % function number = getnumber(interfaces,k)
            
            if nargin==1
                n = numel(interfaces);
                number = cell(0,n);
                for k=1:n
                    number{k} = getnumber(interfaces.interfaces{k});
                end
            else
                number = [];
                for j=1:length(k)
                    number = [number,getnumber(interfaces.interfaces{k(j)})];
                end
            end
        end
        
        function interfaces = setnumber(interfaces,k,l)
            % function interfaces = setnumber(interfaces,k)
            % reset number to interfaces #k
            % throw an error if two interfaces have the same number
            % k is of size numel(interfaces)
            %
            % function interfaces = setnumber(interfaces,k,l)
            % reset number to interfaces #k with #l
            %
            % function interfaces = setnumber(interfaces)
            % set number to interfaces without number
            
            n = numel(interfaces);
            if nargin==1
                number = getnumber(interfaces);
                for k=1:n
                    if isempty(number{k})
                        number{k} = newnumber(interfaces);
                        interfaces.interfaces{k} = setnumber(interfaces.interfaces{k},number{k});
                    end
                end
            elseif nargin==2
                if isa(k,'double')
                    k = num2cell(k);
                end
                if length(k)~=n
                    error('All the interfaces must be renumbered.')
                end
                for j=1:n
                    interfaces.interfaces{j} = setnumber(interfaces.interfaces{j},k(j));
                end 
            elseif nargin==3
                if length(k)~=length(l)
                    error('k and l must have the same size.')
                end
                for j=1:length(k)
                    interfaces.interfaces{k(j)} = setnumber(interfaces.interfaces{k(j)},l(j));
                end  
            else
                error('Wrong bumber of input arguments.')
            end
        end
        
        function number = newnumber(interfaces,excluded)
            % function number = newnumber(interfaces,excluded)
            
            if nargin==1
                excluded = {};
            else
                if ~isa(excluded,'cell')
                    excluded = num2cell(excluded{:});
                end
            end
            allnumber = [getnumber(interfaces),excluded];
            allnumber = [allnumber{:}];
            number = min(setdiff(1:length(allnumber)+1,allnumber));
        end
        
        function checknumber(interfaces)
            % function checknumber(interfaces)
            % throw an error if two interfaces have the same number
            
            number = getnumber(interfaces);
            if length(unique([number{:}]))~=length([number{:}])
                error('The interfaces must have different numbers.')
            end
        end
        
        function [rep,loc] = ismember(k,interfaces)
            % function [rep,loc] = ismember(k,interfaces)
            
            if ~isa(k,'double')
                k = getnumber(k);
                if isa(k,'cell')
                    k = [k{:}];
                end
            end
            if ~isa(interfaces,'double')
                interfaces = getnumber(interfaces);
                if isa(interfaces,'cell')
                    interfaces = [interfaces{:}];
                end
            end
            [rep,loc] = ismember(k,interfaces);
        end
        
        function type = gettype(interfaces,k)
            % function type = gettype(interfaces,k)
            
            if nargin==1
                type = cell(0,interfaces.number);
                for k=1:interfaces.number
                    type{k} = gettype(interfaces.interfaces{k});
                end
            else
                type = [];
                for j=1:length(k)
                   type = [type,gettype(interfaces.interfaces{k(j)})];
                end
            end
        end
        
        function interfaces = horzcat(varargin)
            % function interfaces = horzcat(varargin)
            
            interfaces = Interfaces(varargin{:});
        end
        
        function interfaces = unique(interfaces)
            % function interfaces = unique(interfaces)
            
            number = getnumber(interfaces);
            [~,ia] = unique([number{:}]);
            interfaces.interfaces = interfaces.interfaces(sort(ia));
        end
        
        function interface = getInterface(interfaces,k)
            % function interface = getInterface(interfaces,k)
            
            if nargin==1
                interface = interfaces.interfaces;
            else
                [~,rep] = ismember(k,interfaces);
                if length(k)~=length(rep) || any(rep==0)
                    error('Wrong number')
                end
                interface = interfaces.interfaces(rep);
            end
            if numel(interface)==1
                interface = interface{1};
            elseif isempty(interface)
                interface = [];
            end
        end

        function plotModel(interfaces,varargin)
            % function plotModel(interfaces)
            % Display models interfaces.S of all interfaces
            %
            % function plotModel(interfaces,k)
            % Display models interfaces.S of k selected interfaces
            
            p = ImprovedInputParser;
            addParameter(p,'legend',true,@islogical);
            addParameter(p,'FontSize',16,@isscalar);
            addParameter(p,'LineWidth',0.5,@isscalar);
            addParameter(p,'Interpreter','latex',@ischar);
            if isempty(varargin) || ~isnumeric(varargin{1})
                parse(p,varargin{:})
            else
                parse(p,varargin{2:end})
            end
            
            options = varargin;
            varargin = delcharin({'legend','FontSize','LineWidth','Interpreter'},varargin);
            
            interfaces = getInterface(interfaces,varargin{:});
            
            if numel(interfaces)==1
                plotModel(interfaces,options{:});
            else
                numbers = getnumber([interfaces{:}]);
                figure('Name',['Meshes of interfaces #' num2str([numbers{:}])])
                % set(gcf,'Name',['Meshes of interfaces #' num2str([numbers{:}])])
                hg = cell(1,numel(interfaces));
                leg = cell(1,numel(interfaces));
                for k=1:numel(interfaces)
                    interface = interfaces{k};
                    h = plot(interface.S,'FaceColor',getfacecolor(interface.number+1),'LineWidth',p.Results.LineWidth);
                    hg{k} = hggroup;
                    set(h(:),'Parent',hg{k});
                    leg{k} = ['$\Gamma_{' num2str(interface.number) '}$'];
                end
                if p.Results.legend
                    legend([hg{:}],leg{:},'Location','NorthEastOutside','Interpreter',p.Results.Interpreter)
                end
                set(gca,'FontSize',p.Results.FontSize)
            end
            
        end
        
    end
    
end
