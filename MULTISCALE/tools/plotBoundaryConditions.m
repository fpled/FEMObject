function varargout = plotBoundaryConditions(S,varargin)
% function varargout = plotBoundaryConditions(S)
% Display boundary conditions associated to model S
% S: MODEL

p = ImprovedInputParser;
addParameter(p,'legend',true,@islogical);
addParameter(p,'FontSize',12,@isscalar);
addParameter(p,'LineWidth',0.5,@isscalar);
addParameter(p,'Interpreter','latex',@ischar);
addParameter(p,'FaceColor','w',@(x) isnumeric(x) || ischar(x));
parse(p,varargin{:})

varargin = delcharin({'legend','FontSize','LineWidth','Interpreter','FaceColor'},varargin);

figure('Name','Boundary conditions')
% set(gcf,'Name','Boundary conditions')
clf
h = plot(S,'FaceColor',p.Results.FaceColor,'LineWidth',p.Results.LineWidth,varargin{:});
hg = hggroup;
set(h(:),'Parent',hg);
hold on
[hD,leg] = plotbcond(S);

if p.Results.legend
    % legend(g,'$\Omega$',hD,leg{:},'Interpreter',p.Results.Interpreter)
    legend(hD,leg{:});
else
    legend('off');
end
set(gca,'FontSize',p.Results.FontSize)

if nargout>=1
    varargout{1} = hD;
end
if nargout>=2
    varargout{2} = leg;
end

end
