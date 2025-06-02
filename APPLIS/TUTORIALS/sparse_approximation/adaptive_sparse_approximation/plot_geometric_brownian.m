function plot_geometric_brownian(Xref,Xpc,varargin)
% function geometric_brownian(Xref,Xpc,varargin)
% Display geometric Brownian motion

p = inputParser;
addParameter(p,'legend',true,@islogical);
addParameter(p,'label',true,@islogical);
addParameter(p,'grid',true,@islogical);
addParameter(p,'box',true,@islogical);
addParameter(p,'FontSize',16,@isscalar);
addParameter(p,'LineWidth',1,@isscalar);
addParameter(p,'Interpreter','latex',@ischar);
parse(p,varargin{:})

figure('Name','Evolution of Brownian motion w.r.t. time')
% set(gcf,'Name','Evolution of Brownian motion w.r.t. time')
clf
t = linspace(0,1,length(Xref));
plot(t,[Xref(:),Xpc(:)],'LineWidth',p.Results.LineWidth);
if p.Results.grid
    grid on
end
if p.Results.box
    box on
end
set(gca,'FontSize',p.Results.FontSize)
xlabel('Time')
ylabel('Brownian motion')
if p.Results.legend
    legend({'$B_t^{\mathrm{ref}}$','$B_t^{\mathrm{approx}}$'},'Interpreter',p.Results.Interpreter)
end

end
