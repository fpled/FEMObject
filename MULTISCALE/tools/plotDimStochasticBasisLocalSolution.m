function plotDimStochasticBasisLocalSolution(output,varargin)
% function plotDimStochasticBasisLocalSolution(output,varargin)
% Display the dimension of the stochastic basis of local solution w

p = ImprovedInputParser;
addParameter(p,'legend',true,@islogical);
addParameter(p,'grid',true,@islogical);
addParameter(p,'box',true,@islogical);
addParameter(p,'FontSize',16,@isscalar);
addParameter(p,'LineWidth',1,@isscalar);
addParameter(p,'Interpreter','latex',@ischar);
parse(p,varargin{:})

n = numel(output.dimBasisLocalSolution);
iter = 1:output.iteration;
figure('Name','Evolution of dimension of stochastic basis w.r.t number of iterations')
% set(gcf,'Name','Evolution of dimension of stochastic basis w.r.t number of iterations')
clf
hold on
for k=1:n
    plot(iter,output.dimBasisLocalSolution{k},'-','Color',getfacecolor(k+1),'LineWidth',p.Results.LineWidth)
    if k==1
        leg = {['$P(w^k_{' num2str(k) '})$']};
    else
        leg = [leg, {['$P(w^k_{' num2str(k) '})$']}];
    end
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
ylabel('Dimension of stochastic basis')
if p.Results.legend
    legend(leg{:},'Interpreter',p.Results.Interpreter)
end

end
