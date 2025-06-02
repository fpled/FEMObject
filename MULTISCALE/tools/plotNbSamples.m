function plotNbSamples(output,varargin)
% function plotNbSamples(output,varargin)
% Display number of samples N

p = ImprovedInputParser;
addParameter(p,'legend',true,@islogical);
addParameter(p,'grid',true,@islogical);
addParameter(p,'box',true,@islogical);
addParameter(p,'FontSize',16,@isscalar);
addParameter(p,'LineWidth',1,@isscalar);
addParameter(p,'Interpreter','latex',@ischar);
parse(p,varargin{:})

figure('Name','Evolution of number of samples w.r.t number of iterations')
% set(gcf,'Name','Evolution of number of samples w.r.t number of iterations')
clf
hold on
iter = 1:output.iteration;
n = numel(output.nbSamples);
leg = cell(1,n);
for k=1:n
    plot(iter,output.nbSamples{k}(iter),'-','Color',getfacecolor(k+1),'LineWidth',p.Results.LineWidth)
    leg{k} = ['$N_{' num2str(k) '}$'];
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
ylabel('Number of samples')
if p.Results.legend
    legend(leg{:},'Interpreter',p.Results.Interpreter)
end

end
