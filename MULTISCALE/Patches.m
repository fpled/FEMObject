classdef Patches
    %
    
    %
    %   A MATLAB class for representing a collection of patches
    %
    %   See also PATCH, INTERFACE, INTERFACES
    
    properties
        patches
    end
    
    methods
        
        % Constructor
        function patches = Patches(varargin)
            % Class Patches
            %
            % function patches = Patches(patches)
            % patches: cell containing objects of type Patch
            %
            % function patches = Patches(n)
            % patches: cell of length n containing objects of type Patch
            
            switch nargin
                case 0
                    patches.patches = {};
                case 1
                    if isa(varargin{1},'Patches')
                        patches = varargin{1};
                    elseif isa(varargin{1},'double')
                        n = varargin{1};
                        patches.patches = cell(1,n);
                        for k=1:n
                            patches.patches{k} = Patch();
                        end
                        patches = setnumber(patches);
                    end
                otherwise
                    patches.patches = cell(1,0);
                    for k=1:nargin
                        if isa(varargin{k},'Patch')
                            patches.patches = [patches.patches, {varargin{k}}];
                        elseif isa(varargin{k},'Patches')
                            patches.patches = [patches.patches, varargin{k}.patch];
                        elseif isa(varargin{k},'cell')
                            patches_k = Patches(varargin{k}{:}) ;
                            patches.patches = [patches.patches, patches_k.patches];
                        end
                    end
                    patches = setnumber(patches);
                    patches = unique(patches);
            end
            
        end
        
        function n = numel(patches)
            % function n = numel(patches)
            
            n = numel(patches.patches);
        end
        
        function number = getnumber(patches,k)
            % function number = getnumber(patches,k)
            
            if nargin==1
                n = numel(patches);
                number = cell(0,n);
                for k=1:n
                    number{k} = getnumber(patches.patches{k});
                end
            else
                number = [];
                for j=1:length(k)
                    number = [number,getnumber(patches.patches{k(j)})];
                end
            end
        end
        
        function patches = setnumber(patches,k,l)
            % function patches = setnumber(patches,k)
            % reset number to patches #k
            % throw an error if two patches have the same number
            % k is of size numel(patches)
            %
            % function patches = setnumber(patches,k,l)
            % reset number to patches #k with #l
            %
            % function patches = setnumber(patches)
            % set number to patches without number
            
            n = numel(patches);
            if nargin==1
                number = getnumber(patches);
                for k=1:n
                    if isempty(number{k})
                        number{k} = newnumber(patches);
                        patches.patches{k} = setnumber(patches.patches{k},number{k});
                    end
                end
            elseif nargin==2
                if isa(k,'double')
                    k = num2cell(k);
                end
                if length(k)~=n
                    error('All the patches must be renumbered.')
                end
                for j=1:n
                    patches.patches{j} = setnumber(patches.patches{j},k(j));
                end
            elseif nargin==3
                if length(k)~=length(l)
                    error('k and l must have the same size.')
                end
                for j=1:length(k)
                    patches.patches{k(j)} = setnumber(patches.patches{k(j)},l(j));
                end
            else
                error('Wrong bumber of input arguments.')
            end
        end
        
        function number = newnumber(patches,excluded)
            % function number = newnumber(patches,excluded)
            
            if nargin==1
                excluded = {};
            else
                if ~isa(excluded,'cell')
                    excluded = num2cell(excluded{:});
                end
            end
            allnumber = [getnumber(patches),excluded];
            allnumber = [allnumber{:}];
            number = min(setdiff(1:length(allnumber)+1,allnumber));
        end
        
        function checknumber(patches)
            % function checknumber(patches)
            % throw an error if two patches have the same number
            
            number = getnumber(patches);
            if length(unique([number{:}]))~=length([number{:}])
                error('The patches must have different numbers.')
            end
        end
        
        function [rep,loc] = ismember(k,patches)
            % function [rep,loc] = ismember(k,patches)
            
            if ~isa(k,'double')
                k = getnumber(k);
                if isa(k,'cell')
                    k = [k{:}];
                end
            end
            if ~isa(patches,'double')
                patches = getnumber(patches);
                if isa(patches,'cell')
                    patches = [patches{:}];
                end
            end
            [rep,loc] = ismember(k,patches);
        end
        
        function type = gettype(patches,k)
            % function type = gettype(patches,k)
            
            if nargin==1
                type = cell(0,patches.number);
                for k=1:patches.number
                    type{k} = gettype(patches.patches{k});
                end
            else
                type = [];
                for j=1:length(k)
                    type = [type,gettype(patches.patches{k(j)})];
                end
            end
        end
        
        function patches = horzcat(varargin)
            % function patches = horzcat(varargin)
            
            patches = Patches(varargin{:});
        end
        
        function patches = unique(patches)
            % function patches = unique(patches)
            
            number = getnumber(patches);
            [~,ia] = unique([number{:}]);
            patches.patches = patches.patches(sort(ia));
        end
        
        function patch = getPatch(patches,k)
            % function patch = getPatch(patches,k)
            
            if nargin==1
                patch = patches.patches;
            else
                [~,rep] = ismember(k,patches);
                if length(k)~=length(rep) || any(rep==0)
                    error('Wrong number')
                end
                patch = patches.patches(rep);
            end
            if numel(patch)==1
                patch = patch{1};
            elseif isempty(patch)
                patch = [];
            end
        end
        
        function patches = eval(patches,x,varargin)
            % function patches = eval(patches,x,varargin)
            
            for k=1:numel(patches)
                patches.patches{k} = eval(patches.patches{k},x,varargin{:});
            end
        end
        
        function patches = linearizeOperator(patches,u)
            % function patches = linearizeOperator(patches,u)
            
            for k=1:numel(patches)
                patches.patches{k} = linearize(patches.patches{k},u{k});
            end
        end
        
        function patches = calcOperator(patches,varargin)
            % function patches = calcOperator(patches,varargin)
            
            for k=1:numel(patches)
                patches.patches{k} = calcOperator(patches.patches{k},varargin{:});
            end
        end
        
        function patches = setBoundaryCondition(patches,interfaces)
            % function patches = setBoundaryCondition(patches,interfaces)
            
            for k=1:numel(patches)
                patches.patches{k} = setBoundaryCondition(patches.patches{k},interfaces.interfaces{k});
            end
        end
        
        function plotModel(patches,varargin)
            % function plotModel(patches)
            % Display models patches.S of all patches
            %
            % function plotModel(patches,k)
            % Display models patches.S of k selected patches
            
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
            
            patches = getPatch(patches,varargin{:});
            
            if numel(patches)==1
                plotModel(patches,options{:});
            else
                numbers = getnumber([patches{:}]);
                figure('Name',['Meshes of patches #' num2str([numbers{:}])])
                % set(gcf,'Name',['Meshes of patches #' num2str([numbers{:}])])
                hg = cell(1,numel(patches));
                leg = cell(1,numel(patches));
                for k=1:numel(patches)
                    patch = patches{k};
                    h = plot(patch.S,'FaceColor',getfacecolor(patch.number+1),'LineWidth',p.Results.LineWidth);
                    hg{k} = hggroup;
                    set(h(:),'Parent',hg{k});
                    leg{k} = ['$\Lambda_{' num2str(patch.number) '}$'];
                end
                if p.Results.legend
                    legend([hg{:}],leg{:},'Location','NorthEastOutside','Interpreter',p.Results.Interpreter)
                end
                set(gca,'FontSize',p.Results.FontSize)
            end
            
        end
        
    end
    
end
