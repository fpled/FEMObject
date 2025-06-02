function plotModel(S,varargin)
% function plotModel(S)
% Display model S
% S: MODEL

p = ImprovedInputParser;
addParameter(p,'legend',true,@islogical);
addParameter(p,'node',false,@isscalar);
addParameter(p,'FontSize',16,@isscalar);
addParameter(p,'LineWidth',0.5,@isscalar);
addParameter(p,'Interpreter','latex',@ischar);
addParameter(p,'FaceColor',getfacecolor(1),@(x) isnumeric(x) || ischar(x));
parse(p,varargin{:})

varargin = delcharin({'legend','node','FontSize','LineWidth','Interpreter','FaceColor'},varargin);

figure('Name','Mesh')
% set(gcf,'Name','Mesh')
clf
if ~p.Results.node
    h = plot(S,'FaceColor',p.Results.FaceColor,'LineWidth',p.Results.LineWidth,varargin{:});
else
    h = plot(S,'FaceColor',p.Results.FaceColor,'LineWidth',p.Results.LineWidth,'node',varargin{:});
end
hg = hggroup;
set(h(:),'Parent',hg);
if p.Results.legend
    legend(hg,'$\Omega$','Interpreter',p.Results.Interpreter)
end
set(gca,'FontSize',p.Results.FontSize)

end
