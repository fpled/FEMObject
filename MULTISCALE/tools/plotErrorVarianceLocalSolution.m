function plotErrorVarianceLocalSolution(patches,w,w_ref,varargin)
% function plotErrorVarianceLocalSolution(patches,w,w_ref,varargin)
% Display the relative error in variance of local solution w
% with respect to reference local solution w_ref
% patches: Patches
% w: FunctionalBasisArray of local solution w
% w_ref: FunctionalBasisArray of reference local solution w

p = ImprovedInputParser;
addParameter(p,'colorbar',true,@(x) islogical(x) || ischar(x));
addParameter(p,'colormap','default',@(x) isnumeric(x) || ischar(x));
addParameter(p,'FontSize',16,@isscalar);
parse(p,varargin{:})

varargin = delcharin({'colorbar','colormap','FontSize'},varargin);
n = numel(patches);

if ischarin('displ',varargin)
    i = getcharin('displ',varargin);
    figure('Name',['Relative error in variance of w_' num2str(i) ' w.r.t. w_ref_' num2str(i)])
    % set(gcf,'Name',['Relative error in variance of w_' num2str(i) ' w.r.t. w_ref_' num2str(i)])
else
    figure('Name','Relative error in variance of w w.r.t. w_ref')
    % set(gcf,'Name','Relative error in variance of n w.r.t. w_ref')
end
clf

for k=1:n
    patch = patches.patches{k};
    err_patch = abs(variance(w{patch.number})-variance(w_ref{patch.number}))'./max(variance(w{patch.number}));
    plot_sol(patch.S,err_patch,varargin{:});
    colormap(p.Results.colormap)
    if p.Results.colorbar
        colorbar
    end
    set(gca,'FontSize',p.Results.FontSize)
end

end

