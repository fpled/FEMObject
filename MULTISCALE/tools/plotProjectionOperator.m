function plotProjectionOperator(glob,patch,numnode,varargin)
% function plotProjectionOperator(glob,patch,numnode)
% Display coupling nodes for projection operator
% glob: Global or GlobalOutside
% patches: Patches
% numnode: double containing the numbers of coupling nodes in global domain if glob is an instance of class Global
%          or in complementary subdomain if glob is an instance of class GlobalOutside

p = ImprovedInputParser;
addParameter(p,'legend',true,@islogical);
addParameter(p,'FontSize',16,@isscalar);
addParameter(p,'LineWidth',1,@isscalar);
addParameter(p,'Interpreter','latex',@ischar);
addParameter(p,'EdgeColor','k',@(x) isnumeric(x) || ischar(x));
parse(p,varargin{:})

if isa(glob,'Global')
    S_out = glob.S_out;
elseif isa(glob,'GlobalOutside')
    S_out = glob.S;
end
    
figure('Name',['Coupling nodes for projection operator over interface #' num2str(patch.number)])
% set(gcf,'Name',['Coupling nodes for projection operator over interface #' num2str(patch.number)])
clf
h_out = plot(S_out,'FaceColor',getfacecolor(1),'EdgeColor',p.Results.EdgeColor,'LineWidth',p.Results.LineWidth);
h_patch = plot(patch.S,'FaceColor',getfacecolor(patch.number+1),'EdgeColor',p.Results.EdgeColor,'LineWidth',p.Results.LineWidth);
h = plot(glob.S.node(numnode),'o','Color','b','MarkerFaceColor','b','LineWidth',p.Results.LineWidth);
hg_out = hggroup;
hg_patch = hggroup;
hg = hggroup;
set(h_out(:),'Parent',hg_out);
set(h_patch(:),'Parent',hg_patch);
set(h(:),'Parent',hg);
if p.Results.legend
    legend([hg_out,hg_patch,h],'$\Omega \setminus \Lambda$',...
        ['$\Lambda_{' num2str(patch.number) '}$'],...
        ['$i \in \partial(\Omega \setminus \Lambda) \cap \partial \Lambda_{' num2str(patch.number) '}$'],...
        'Location','NorthEastOutside','Interpreter',p.Results.Interpreter)
end
set(gca,'FontSize',p.Results.FontSize)

end
