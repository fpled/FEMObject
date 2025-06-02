function plotDimStochasticBasis(output,varargin)
% function plotDimStochasticBasis(output,varargin)
% Display the dimension of the stochastic basis of local solution w and
% Lagrange Multiplier lambda 

p = ImprovedInputParser;
addParameter(p,'legend',true,@islogical);
addParameter(p,'grid',true,@islogical);
addParameter(p,'box',true,@islogical);
addParameter(p,'FontSize',16,@isscalar);
addParameter(p,'LineWidth',1,@isscalar);
addParameter(p,'Interpreter','latex',@ischar);
parse(p,varargin{:})

figure('Name','Evolution of dimension of stochastic basis w.r.t number of iterations')
% set(gcf,'Name','Evolution of dimension of stochastic basis w.r.t number of iterations')
clf
hold on
iter = 1:output.iteration;
leg = {};
n = numel(output.dimBasisLocalSolution);
for k=1:n
    plot(iter,output.dimBasisLocalSolution{k},'-','Color',getfacecolor(k+1),'LineWidth',p.Results.LineWidth)
    leg = [leg, {['$P(w^k_{' num2str(k) '})$']}];
end
n = numel(output.dimBasisLagrangeMultiplier);
for k=1:n
    plot(iter,output.dimBasisLagrangeMultiplier{k},'--','Color',getfacecolor(k+1),'LineWidth',p.Results.LineWidth)
    leg = [leg, {['$P(\lambda^k_{' num2str(k) '})$']}];
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
