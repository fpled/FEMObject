function plotStagnationGlobalSolution(output,varargin)
% function plotStagnationGlobalSolution(output,varargin)
% Display the stagnation indicator of global solution U

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
