function plotModelDeflection(S,u,varargin)
% function plotModelDeflection(S,u,ampl)
% Display deflected model S
% S: MODEL
% u: solution
% ampl: amplification factor

p = ImprovedInputParser;
addParameter(p,'legend',true,@islogical);
addParameter(p,'FontSize',16,@isscalar);
addParameter(p,'LineWidth',0.5,@isscalar);
addParameter(p,'ampl',getsize(S)/max(abs(unfreevector(S,u))),@isscalar);
addParameter(p,'Interpreter','latex',@ischar);
addParameter(p,'FaceColor',getfacecolor(1),@(x) isnumeric(x) || ischar(x));
parse(p,varargin{:})

varargin = delcharin({'legend','FontSize','LineWidth','Interpreter','FaceColor','ampl'},varargin);

u = unfreevector(S,u);

figure('Name','Deflection of mesh')
% set(gcf,'Name','Deflection of mesh')
clf
h = plot(S+p.Results.ampl*u,'FaceColor',p.Results.FaceColor,'LineWidth',p.Results.LineWidth,varargin{:});
hg = hggroup;
set(h(:),'Parent',hg);
if p.Results.legend
    legend(hg,'$\Omega$','Interpreter',p.Results.Interpreter)
end
set(gca,'FontSize',p.Results.FontSize)

end
