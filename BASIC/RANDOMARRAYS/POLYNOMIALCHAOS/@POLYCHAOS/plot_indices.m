function plot_indices(PC,varargin)
% function plot_indices(PC,'dim',dim)
% Plot the multi-index set of POLYCHAOS PC projected on dimension dim
% dim = 1:min(3,PC.M) by default
% 
% function plot_indices(PC,'dim',dim,'MarkerArea',MarkerArea,'MarkerColor',MarkerColor,'MarkerType',MarkerType,'axis',axis)
% Plot the multi-index set of POLYCHAOS PC projected on dimension dim with
% optional parameters
% marker = 'o' by default
% color = 'b' by default
% 
% function plot_indices(PC,'dim',dim,'maximal',true)
% Plot the multi-index set of POLYCHAOS PC projected on dimension dim
% and its maximal indices
% 
% function plot_indices(PC,'dim',dim,'margin',true)
% Plot the multi-index set of POLYCHAOS PC projected on dimension dim
% and its margin
% 
% function plot_indices(PC,'dim',dim,'reducedmargin',true)
% Plot the multi-index set of POLYCHAOS PC projected on dimension dim
% and its reduced margin

isboolean = @(x) islogical(x) || isnumeric(x) && all(x(:)==0 | x(:)==1);

p = inputParser;
addParameter(p,'dim',1:min(3,PC.M),@isnumeric);
addParameter(p,'MarkerArea',36,@isnumeric);
addParameter(p,'MarkerColor','b',@(x) isnumeric(x) || ischar(x));
addParameter(p,'MarkerType','o',@ischar);
addParameter(p,'MarkerArea_max',216,@isnumeric);
addParameter(p,'MarkerType_max','s',@ischar);
addParameter(p,'MarkerColor_marg','r',@(x) isnumeric(x) || ischar(x));
addParameter(p,'legend',true,isboolean);
addParameter(p,'label',true,isboolean);
addParameter(p,'grid',true,isboolean);
addParameter(p,'axis',[],@isnumeric)
addParameter(p,'tick',[],@(x) isnumeric(x) || iscell(x))
addParameter(p,'maximal',false,isboolean);
addParameter(p,'margin',false,isboolean);
addParameter(p,'reducedMargin',false,isboolean);
addParameter(p,'FontSize',16,@isscalar);
addParameter(p,'LineWidth',1,@isscalar);
addParameter(p,'Interpreter','latex',@ischar);
parse(p,varargin{:})

PC_dim = projectdim(PC,p.Results.dim);
ind = getindices(PC);
ind = unique(ind(:,p.Results.dim),'rows');

switch length(p.Results.dim)
    case 0
        error('The number of dimensions must be greater or equal to 1!')
    case 1
        scatter(ind(:,1),zeros(size(ind,1),1),p.Results.MarkerArea,p.Results.MarkerColor,p.Results.MarkerType,'LineWidth',p.Results.LineWidth,'MarkerEdgeColor',p.Results.MarkerColor,'MarkerFaceColor',p.Results.MarkerColor)
        axis equal
        set(gca,'FontSize',p.Results.FontSize)
        if ~isempty(p.Results.tick)
            set(gca,'XTick',p.Results.tick,'YTick',[])
        else
            set(gca,'XTick',0:max(ind(:,1)),'YTick',[])
        end
        if ~isempty(p.Results.axis)
            axis(p.Results.axis)
        else
            set(gca,'XLim',[-eps,max(ind(:,1))+1],'YLim',[-eps,1])
        end
        leg = {['$\alpha \in \mathcal{A}_{' num2str(length(PC)) '}$']};
        hold on
        
        if p.Results.maximal
            ind_max = getindices(PC_dim,'max');
            scatter(ind_max(:,1),zeros(size(ind_max,1),1),p.Results.MarkerArea_max,p.Results.MarkerColor,p.Results.MarkerType_max,'LineWidth',p.Results.LineWidth)
            leg = [leg, {['$\alpha \in \mathcal{M}(\mathcal{A}_{' num2str(length(PC)) '})$']}];
        end
        
        if p.Results.margin
            ind_marg = getindices(PC_dim,'margin');
            scatter(ind_marg(:,1),zeros(size(ind_marg,1),1),p.Results.MarkerArea,p.Results.MarkerColor_marg,p.Results.MarkerType,'LineWidth',p.Results.LineWidth,'MarkerEdgeColor',p.Results.MarkerColor_marg,'MarkerFaceColor',p.Results.MarkerColor_marg)
            if ~isempty(p.Results.tick)
                set(gca,'XTick',p.Results.tick,'YTick',[])
            else
                set(gca,'XTick',0:max(ind_marg(:,1)),'YTick',[])
            end
            if ~isempty(p.Results.axis)
                axis(p.Results.axis)
            else
                set(gca,'XLim',[-eps,max(ind_marg(:,1))+1],'YLim',[-eps,1])
            end
            leg = [leg, {['$\alpha \in \mathcal{N}(\mathcal{A}_{' num2str(length(PC)) '})$']}];
            if p.Results.reducedMargin
                ind_red_marg = getindices(PC_dim,'reducedMargin');
                scatter(ind_red_marg(:,1),zeros(size(ind_red_marg,1),1),p.Results.MarkerArea_max,p.Results.MarkerColor_marg,p.Results.MarkerType_max,'LineWidth',p.Results.LineWidth)
                leg = [leg, {['$\alpha \in \mathcal{N}_r(\mathcal{A}_{' num2str(length(PC)) '})$']}];
            end
        elseif p.Results.reducedMargin
            ind_red_marg = getindices(PC_dim,'reducedMargin');
            scatter(ind_red_marg(:,1),zeros(size(ind_red_marg,1),1),p.Results.MarkerArea,p.Results.MarkerColor_marg,p.Results.MarkerType,'LineWidth',p.Results.LineWidth,'MarkerEdgeColor',p.Results.MarkerColor_marg,'MarkerFaceColor',p.Results.MarkerColor_marg)
            if ~isempty(p.Results.tick)
                set(gca,'XTick',p.Results.tick,'YTick',[])
            else
                set(gca,'XTick',0:max(ind_red_marg(:,1)),'YTick',[])
            end
            if ~isempty(p.Results.axis)
                axis(p.Results.axis)
            else
                set(gca,'XLim',[-eps,max(ind_red_marg(:,1))+1],'YLim',[-eps,1])
            end
            leg = [leg, {['$\alpha \in \mathcal{N}_r(\mathcal{A}_{' num2str(length(PC)) '})$']}];
        end
        if p.Results.label
            xlabel(['$\alpha_{' num2str(p.Results.dim) '}$'],'Interpreter',p.Results.Interpreter)
        end
        hold off
    case 2
        scatter(ind(:,1),ind(:,2),p.Results.MarkerArea,p.Results.MarkerColor,p.Results.MarkerType,'LineWidth',p.Results.LineWidth,'MarkerEdgeColor',p.Results.MarkerColor,'MarkerFaceColor',p.Results.MarkerColor)
        axis equal
        set(gca,'FontSize',p.Results.FontSize)
        if ~isempty(p.Results.tick)
            set(gca,'XTick',p.Results.tick{1},'YTick',p.Results.tick{2})
        else
            set(gca,'XTick',0:max(ind(:,1)),'YTick',0:max(ind(:,2)))
        end
        if ~isempty(p.Results.axis)
            axis(p.Results.axis)
        else
            set(gca,'XLim',[-eps,max(ind(:,1))+1],'YLim',[-eps,max(ind(:,2))+1])
        end
        leg = {['$\alpha \in \mathcal{A}_{' num2str(length(PC)) '}$']};
        hold on
        if p.Results.maximal
            ind_max = getindices(PC_dim,'max');
            scatter(ind_max(:,1),ind_max(:,2),p.Results.MarkerArea_max,p.Results.MarkerColor,p.Results.MarkerType_max,'LineWidth',p.Results.LineWidth)
            leg = [leg, {['$\alpha \in \mathcal{M}(\mathcal{A}_{' num2str(length(PC)) '})$']}];
        end
        if p.Results.margin
            ind_marg = getindices(PC_dim,'margin');
            scatter(ind_marg(:,1),ind_marg(:,2),p.Results.MarkerArea,p.Results.MarkerColor_marg,p.Results.MarkerType,'LineWidth',p.Results.LineWidth,'MarkerEdgeColor',p.Results.MarkerColor_marg,'MarkerFaceColor',p.Results.MarkerColor_marg)
            if ~isempty(p.Results.tick)
                set(gca,'XTick',p.Results.tick{1},'YTick',p.Results.tick{2})
            else
                set(gca,'XTick',0:max(ind_marg(:,1)),'YTick',0:max(ind_marg(:,2)))
            end
            if ~isempty(p.Results.axis)
                axis(p.Results.axis)
            else
                set(gca,'XLim',[-eps,max(ind_marg(:,1))+1],'YLim',[-eps,max(ind_marg(:,2))+1])
            end
            leg = [leg, {['$\alpha \in \mathcal{N}(\mathcal{A}_{' num2str(length(PC)) '})$']}];
            if p.Results.reducedMargin
                ind_red_marg = getindices(PC_dim,'reducedMargin');
                scatter(ind_red_marg(:,1),ind_red_marg(:,2),p.Results.MarkerArea_max,p.Results.MarkerColor_marg,p.Results.MarkerType_max,'LineWidth',p.Results.LineWidth)
                leg = [leg, {['$\alpha \in \mathcal{N}_r(\mathcal{A}_{' num2str(length(PC)) '})$']}];
            end
        elseif p.Results.reducedMargin
            ind_red_marg = getindices(PC_dim,'reducedMargin');
            scatter(ind_red_marg(:,1),ind_red_marg(:,2),p.Results.MarkerArea,p.Results.MarkerColor_marg,p.Results.MarkerType,'LineWidth',p.Results.LineWidth,'MarkerEdgeColor',p.Results.MarkerColor_marg,'MarkerFaceColor',p.Results.MarkerColor_marg)
            if ~isempty(p.Results.tick)
                set(gca,'XTick',p.Results.tick{1},'YTick',p.Results.tick{2})
            else
                set(gca,'XTick',0:max(ind_red_marg(:,1)),'YTick',0:max(ind_red_marg(:,2)))
            end
            if ~isempty(p.Results.axis)
                axis(p.Results.axis)
            else
                set(gca,'XLim',[-eps,max(ind_red_marg(:,1))+1],'YLim',[-eps,max(ind_red_marg(:,2))+1])
            end
            leg = [leg, {['$\alpha \in \mathcal{N}_r(\mathcal{A}_{' num2str(length(PC)) '})$']}];
        end
        if p.Results.label
            xlabel(['$\alpha_{' num2str(p.Results.dim(1)) '}$'],'Interpreter',p.Results.Interpreter)
            ylabel(['$\alpha_{' num2str(p.Results.dim(2)) '}$'],'Interpreter',p.Results.Interpreter)
        end
        hold off
    case 3
        scatter3(ind(:,1),ind(:,2),ind(:,3),p.Results.MarkerArea,p.Results.MarkerColor,p.Results.MarkerType,'LineWidth',p.Results.LineWidth,'MarkerEdgeColor',p.Results.MarkerColor,'MarkerFaceColor',p.Results.MarkerColor)
        axis equal
        % view(37.5,30)
        view(50,10)
        set(gca,'FontSize',p.Results.FontSize)
        if ~isempty(p.Results.tick)
            set(gca,'XTick',p.Results.tick{1},'YTick',p.Results.tick{2},'ZTick',p.Results.tick{3})
        else
            set(gca,'XTick',0:max(ind(:,1)),'YTick',0:max(ind(:,2)),'ZTick',0:max(ind(:,3)))
        end
        if ~isempty(p.Results.axis)
            axis(p.Results.axis)
        else
            set(gca,'XLim',[-eps,max(ind(:,1))+1],'YLim',[-eps,max(ind(:,2))+1],'ZLim',[-eps,max(ind(:,3))+1])
        end
        leg = {['$\alpha \in \mathcal{A}_{' num2str(length(PC)) '}$']};
        hold on
        if p.Results.maximal
            ind_max = getindices(PC_dim,'max');
            scatter3(ind_max(:,1),ind_max(:,2),ind_max(:,3),p.Results.MarkerArea*4,p.Results.MarkerColor,p.Results.MarkerType_max,'LineWidth',p.Results.LineWidth)
            leg = [leg, {['$\alpha \in \mathcal{M}(\mathcal{A}_{' num2str(length(PC)) '})$']}];
        end
        if p.Results.margin
            ind_marg = getindices(PC_dim,'margin');
            scatter3(ind_marg(:,1),ind_marg(:,2),ind_marg(:,3),p.Results.MarkerArea,p.Results.MarkerColor_marg,p.Results.MarkerType,'LineWidth',p.Results.LineWidth,'MarkerEdgeColor',p.Results.MarkerColor_marg,'MarkerFaceColor',p.Results.MarkerColor_marg)
            if ~isempty(p.Results.tick)
                set(gca,'XTick',p.Results.tick{1},'YTick',p.Results.tick{2},'ZTick',p.Results.tick{3})
            else
                set(gca,'XTick',0:max(ind_marg(:,1)),'YTick',0:max(ind_marg(:,2)),'ZTick',0:max(ind_marg(:,3)))
            end
            if ~isempty(p.Results.axis)
                axis(p.Results.axis)
            else
                set(gca,'XLim',[-eps,max(ind_marg(:,1))+1],'YLim',[-eps,max(ind_marg(:,2))+1],'ZLim',[-eps,max(ind_marg(:,3))+1])
            end
            leg = [leg, {['$\alpha \in \mathcal{N}(\mathcal{A}_{' num2str(length(PC)) '})$']}];
            if p.Results.reducedMargin
                ind_red_marg = getindices(PC_dim,'reducedMargin');
                scatter3(ind_red_marg(:,1),ind_red_marg(:,2),ind_red_marg(:,3),p.Results.MarkerArea_max,p.Results.MarkerColor_marg,p.Results.MarkerType_max,'LineWidth',p.Results.LineWidth)
                leg = [leg, {['$\alpha \in \mathcal{N}_r(\mathcal{A}_{' num2str(length(PC)) '})$']}];
            end
        elseif p.Results.reducedMargin
            ind_red_marg = getindices(PC_dim,'reducedMargin');
            scatter3(ind_red_marg(:,1),ind_red_marg(:,2),ind_red_marg(:,3),p.Results.MarkerArea,p.Results.MarkerColor_marg,p.Results.MarkerType,'LineWidth',p.Results.LineWidth,'MarkerEdgeColor',p.Results.MarkerColor_marg,'MarkerFaceColor',p.Results.MarkerColor_marg)
            if ~isempty(p.Results.tick)
                set(gca,'XTick',p.Results.tick{1},'YTick',p.Results.tick{2},'ZTick',p.Results.tick{3})
            else
                set(gca,'XTick',0:max(ind_red_marg(:,1)),'YTick',0:max(ind_red_marg(:,2)),'ZTick',0:max(ind_red_marg(:,3)))
            end
            if ~isempty(p.Results.axis)
                axis(p.Results.axis)
            else
                set(gca,'XLim',[-eps,max(ind_red_marg(:,1))+1],'YLim',[-eps,max(ind_red_marg(:,2))+1],'ZLim',[-eps,max(ind_red_marg(:,3))+1])
            end
            leg = [leg, {['$\alpha \in \mathcal{N}_r(\mathcal{A}_{' num2str(length(PC)) '})$']}];
        end
        if p.Results.label
            xlabel(['$\alpha_{' num2str(p.Results.dim(1)) '}$'],'Interpreter',p.Results.Interpreter)
            ylabel(['$\alpha_{' num2str(p.Results.dim(2)) '}$'],'Interpreter',p.Results.Interpreter)
            zlabel(['$\alpha_{' num2str(p.Results.dim(3)) '}$'],'Interpreter',p.Results.Interpreter)
        end
        hold off
    otherwise
        error('The number of dimensions must not exceed 3 for display!')
end
if p.Results.grid
    grid on
end
if p.Results.legend
    legend(leg{:},'Interpreter',p.Results.Interpreter)
end

end
