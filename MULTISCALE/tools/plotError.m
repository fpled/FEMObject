function plotError(output,varargin)
% function plotError(output,varargin)
% Display the error indicator of global solution U, local solution w and
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

figure('Name','Evolution of error indicators w.r.t number of iterations')
% set(gcf,'Name','Evolution of error indicators w.r.t number of iterations')
clf
iter = 1:output.iteration;
semilogy([0,iter],[output.errorGlobalSolutionInit,output.errorGlobalSolution(iter)],'-k','LineWidth',p.Results.LineWidth)
if israndom
    leg = {'$\varepsilon_{\Xi,\Omega\setminus\Lambda}(U^k;U^{\mathrm{ref}})$'};
else
    leg = {'$\varepsilon_{\Omega\setminus\Lambda}(U^k;U^{\mathrm{ref}})$'};
end
hold on
n = numel(output.errorLocalSolution);
for k=1:n
    semilogy([0,iter],[output.errorLocalSolutionInit{k},output.errorLocalSolution{k}(iter)],'-','Color',getfacecolor(k+1),'LineWidth',p.Results.LineWidth)
    if israndom
        leg = [leg, {['$\varepsilon_{\Xi,\Lambda_{' num2str(k) '}}(w^k_{' num2str(k) '};w^{\mathrm{ref}}_{' num2str(k) '})$']}];
    else
        leg = [leg, {['$\varepsilon_{\Lambda_{' num2str(k) '}}(w^k_{' num2str(k) '};w^{\mathrm{ref}}_{' num2str(k) '})$']}];
    end
end
% n = numel(output.errorLagrangeMultiplier);
for k=1:n
    semilogy([0,iter],[output.errorLagrangeMultiplierInit{k},output.errorLagrangeMultiplier{k}(iter)],'--','Color',getfacecolor(k+1),'LineWidth',p.Results.LineWidth)
    if israndom
        leg = [leg, {['$\varepsilon_{\Xi,\Gamma_{' num2str(k) '}}(\lambda^k_{' num2str(k) '};\lambda^{\mathrm{ref}}_{' num2str(k) '})$']}];
    else
        leg = [leg, {['$\varepsilon_{\Gamma_{' num2str(k) '}}(\lambda^k_{' num2str(k) '};\lambda^{\mathrm{ref}}_{' num2str(k) '})$']}];
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
ylabel('Error indicator')
if p.Results.legend
    legend(leg{:},'Interpreter',p.Results.Interpreter)
end

end
