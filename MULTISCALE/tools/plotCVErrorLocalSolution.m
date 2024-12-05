function plotCVErrorLocalSolution(output,varargin)
% function plotCVErrorLocalSolution(output,varargin)
% Display the cross-validation error of local solution w

p = ImprovedInputParser;
addParameter(p,'legend',true,@islogical);
addParameter(p,'grid',true,@islogical);
addParameter(p,'box',true,@islogical);
addParameter(p,'FontSize',16,@isscalar);
addParameter(p,'LineWidth',1,@isscalar);
addParameter(p,'Interpreter','latex',@ischar);
parse(p,varargin{:})

figure('Name','Evolution of cross-validation error w.r.t number of iterations')
% set(gcf,'Name','Evolution of cross-validation error w.r.t number of iterations')
clf
hold on
iter = 1:output.iteration;
leg = {};
n = numel(output.CVErrorLocalSolution);
for k=1:n
    semilogy(iter,cell2mat(cellfun(@norm,output.CVErrorLocalSolution{k},'UniformOutput',false)),'-','Color',getfacecolor(k+1),'LineWidth',p.Results.LineWidth)
    leg = [leg, {['$\varepsilon_{\mathrm{cv}}(w^k_{' num2str(k) '})$']}];
end
hold off
if p.Results.grid
    grid on
end
if p.Results.box
    box on
end
set(gca,'FontSize',p.Results.FontSize)
xlabel('Number of iterations')
ylabel('Cross-validation error')
if p.Results.legend
    l = legend(leg{:});
    set(l,'Interpreter',p.Results.Interpreter)
end

end
