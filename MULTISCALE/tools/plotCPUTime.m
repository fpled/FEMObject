function plotCPUTime(output,varargin)
% function plotCPUTime(output,varargin)
% Display CPU time

p = ImprovedInputParser;
addParameter(p,'legend',true,@islogical);
addParameter(p,'grid',true,@islogical);
addParameter(p,'box',true,@islogical);
addParameter(p,'FontSize',16,@isscalar);
addParameter(p,'LineWidth',1,@isscalar);
parse(p,varargin{:})

figure('Name','Evolution of CPU time w.r.t. number of iterations')
% set(gcf,'Name','Evolution of CPU time w.r.t. number of iterations')
clf
iter = 1:output.iteration;
plot(iter,output.time(iter),'-r','LineWidth',p.Results.LineWidth)
if p.Results.grid
    grid on
end
if p.Results.box
    box on
end
set(gca,'FontSize',p.Results.FontSize)
xlabel('Number of iterations')
ylabel('CPU time (s)')
if p.Results.legend
    legend('CPU time versus Nb iterations')
end

end
