function plotErrorGlobalSolution(output,varargin)
% function plotErrorGlobalSolution(output,varargin)
% Display the error indicator of global solution U

p = ImprovedInputParser;
addParameter(p,'legend',true,@islogical);
addParameter(p,'grid',true,@islogical);
addParameter(p,'box',true,@islogical);
addParameter(p,'FontSize',16,@isscalar);
addParameter(p,'LineWidth',1,@isscalar);
addParameter(p,'Interpreter','latex',@ischar);
parse(p,varargin{:})

israndom = isfield(output,'nbSamples');

figure('Name','Evolution of error indicator w.r.t number of iterations')
% set(gcf,'Name','Evolution of error indicator w.r.t number of iterations')
clf
iter = 1:output.iteration;
semilogy([0,iter],[output.errorGlobalSolutionInit,output.errorGlobalSolution(iter)],'-k','LineWidth',p.Results.LineWidth)
if israndom
    leg = {'$\varepsilon_{\Xi,\Omega\setminus\Lambda}(U^k;U^{\mathrm{ref}})$'};
else
    leg = {'$\varepsilon_{\Omega\setminus\Lambda}(U^k;U^{\mathrm{ref}})$'};
end
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
