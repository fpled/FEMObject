function plotRelaxationParameter(output,varargin)
% function plotRelaxationParameter(output,varargin)
% Display relaxation parameter rho

p = ImprovedInputParser;
addParameter(p,'legend',true,@islogical);
addParameter(p,'grid',true,@islogical);
addParameter(p,'box',true,@islogical);
addParameter(p,'FontSize',16,@isscalar);
addParameter(p,'LineWidth',1,@isscalar);
parse(p,varargin{:})

figure('Name','Evolution of relaxation parameter w.r.t. number of iterations')
% set(gcf,'Name','Evolution of relaxation parameter w.r.t. number of iterations')
clf
iter = 1:output.iteration;
plot(iter,output.relaxationParameter(iter),'-b','LineWidth',p.Results.LineWidth)
if p.Results.grid
    grid on
end
if p.Results.box
    box on
end
set(gca,'FontSize',p.Results.FontSize)
xlabel('Number of iterations')
ylabel('Relaxation parameter')
if p.Results.legend
    legend('Relaxation parameter versus Nb iterations')
end

end
