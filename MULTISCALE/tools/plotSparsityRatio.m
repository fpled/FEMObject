function plotSparsityRatio(output,varargin)
% function plotSparsityRatio(output,varargin)
% Display sparsity ratio

p = ImprovedInputParser;
addParameter(p,'legend',true,@islogical);
addParameter(p,'grid',true,@islogical);
addParameter(p,'box',true,@islogical);
addParameter(p,'FontSize',16,@isscalar);
addParameter(p,'LineWidth',1,@isscalar);
addParameter(p,'Interpreter','latex',@ischar);
parse(p,varargin{:})

figure('Name','Evolution of sparsity ratio w.r.t number of iterations')
% set(gcf,'Name','Evolution of sparsity ratio w.r.t number of iterations')
clf
iter = 1:output.iteration;
plot(iter,output.sparsityRatio_U(iter),'-k','LineWidth',p.Results.LineWidth)
leg = {'$\iota(U^k)$'};
hold on
n = numel(output.sparsityRatio_w);
for k=1:n
    semilogy(iter,output.sparsityRatio_w{k}(iter),'-','Color',getfacecolor(k+1),'LineWidth',p.Results.LineWidths)
    leg = [leg, {['$\iota(w^k_{' num2str(k) '})$']}];
end
% n = numel(output.sparsityRatio_lambda);
for k=1:n
    semilogy(iter,output.sparsityRatio_lambda{k}(iter),'--','Color',getfacecolor(k+1),'LineWidth',p.Results.LineWidth)
    leg = [leg, {['$\iota(\lambda^k_{' num2str(k) '})$']}];
end
hold off
if p.Results.grid
    grid on
end
if p.Results.box
    box on
end
set(gca,'FontSize',p.Results.FontSize)
ylim([0 1.1]);
xlabel('Number of iterations')
ylabel('Sparsity ratio')
if p.Results.legend
    legend(leg{:},'Interpreter',p.Results.Interpreter)
end

end
