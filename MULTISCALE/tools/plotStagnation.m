function plotStagnation(output,varargin)
% function plotStagnation(output,varargin)
% Display the stagnation indicator of global solution U, local solution w and
% Lagrange multiplier lambda

p = ImprovedInputParser;
addParameter(p,'legend',true,@islogical);
addParameter(p,'grid',true,@islogical);
addParameter(p,'box',true,@islogical);
addParameter(p,'FontSize',16,@isscalar);
addParameter(p,'LineWidth',1,@isscalar);
addParameter(p,'Interpreter','latex',@ischar);
parse(p,varargin{:})

israndom = isfield(output,'nbSamples');

figure('Name','Evolution of stagnation indicators w.r.t number of iterations')
% set(gcf,'Name','Evolution of stagnation indicators w.r.t number of iterations')
clf
iter = 1:output.iteration;
semilogy(iter,output.stagnationGlobalSolution(iter),'-k','LineWidth',p.Results.LineWidth)
if israndom
    leg = {'$\varepsilon_{\Xi,\Omega}(U^k;U^{k-1})$'};
else
    leg = {'$\varepsilon_{\Omega}(U^k;U^{k-1})$'};
end
hold on
n = numel(output.stagnationLocalSolution);
for k=1:n
    semilogy(iter,output.stagnationLocalSolution{k}(iter),'-','Color',getfacecolor(k+1),'LineWidth',p.Results.LineWidth)
    if israndom
        leg = [leg, {['$\varepsilon_{\Xi,\Lambda_{' num2str(k) '}}(w^k_{' num2str(k) '};w^{k-1}_{' num2str(k) '})$']}];
    else
        leg = [leg, {['$\varepsilon_{\Lambda_{' num2str(k) '}}(w^k_{' num2str(k) '};w^{k-1}_{' num2str(k) '})$']}];
    end
end
% n = numel(output.stagnationLagrangeMultiplier);
for k=1:n
    semilogy(iter,output.stagnationLagrangeMultiplier{k}(iter),'--','Color',getfacecolor(k+1),'LineWidth',p.Results.LineWidth)
    if israndom
        leg = [leg, {['$\varepsilon_{\Xi,\Gamma_{' num2str(k) '}}(\lambda^k_{' num2str(k) '};\lambda^{k-1}_{' num2str(k) '})$']}];
    else
        leg = [leg, {['$\varepsilon_{\Gamma_{' num2str(k) '}}(\lambda^k_{' num2str(k) '};\lambda^{k-1}_{' num2str(k) '})$']}];
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
ylabel('Stagnation indicator')
if p.Results.legend
    legend(leg{:},'Interpreter',p.Results.Interpreter)
end

end
