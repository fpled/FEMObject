function plotSensitivityIndicesLocalSolution(patches,w,alpha,varargin)
% function plotSensitivityIndicesLocalSolution(patches,w,alpha,varargin)
% Display the sensitivity indices of local solution w associated with the
% group of variables alpha in {1,..,d} with d = ndims(w)
% patches: Patches or Patch
% w: FunctionalBasisArray of local solution w
% alpha: 1-by-s array of integers or 1-by-d logical
% - if alpha is an array of integers, indices with respect
% to variables alpha
% - if alpha is logical, indices with respect
% to variables find(alpha)

p = ImprovedInputParser;
addParameter(p,'colorbar',true,@(x) islogical(x) || ischar(x));
addParameter(p,'colormap','default',@(x) isnumeric(x) || ischar(x));
addParameter(p,'FontSize',16,@isscalar);
parse(p,varargin{:})

varargin = delcharin({'colorbar','colormap','FontSize'},varargin);

if isa(patches,'Patches')
    numbers = getnumber(patches);
    if ischarin('displ',varargin)
        i = getcharin('displ',varargin);
        figure('Name',['Sensitivity index of w_' num2str(i) ' over patches #' num2str([numbers{:}]) ' for random variables #' num2str(alpha)])
        % set(gcf,'Name',['Sensitivity index of w_' num2str(i) ' over patches #' num2str([numbers{:}]) ' for random variables #' num2str(alpha)])
    else
        figure('Name',['Sensitivity index of w over patches #' num2str([numbers{:}]) ' for random variables #' num2str(alpha)])
        % set(gcf,'Name',['Sensitivity index of w over patches #' num2str([numbers{:}] ' for random variables #' num2str(alpha)])
    end
    clf
    n = numel(patches);
    for k=1:n
        patch = patches.patches{k};
        v = variance(w{patch.number});
        s = varianceConditionalExpectation(w{patch.number},alpha)'./max(v);
        plot_sol(patch.S,s,varargin{:});
        colormap(p.Results.colormap)
        if p.Results.colorbar
            colorbar
        end
        set(gca,'FontSize',p.Results.FontSize)
    end
elseif isa(patches,'Patch')
    patch = patches;
    if ischarin('displ',varargin)
        i = getcharin('displ',varargin);
        figure('Name',['Sensitivity index of w_' num2str(i) ' over patch #' num2str(patch.number) ' for random variables #' num2str(alpha)])
        % set(gcf,'Name',['Sensitivity index of w_' num2str(i) ' over patch #' num2str(patch.number) ' for random variables #' num2str(alpha)])
    else
        figure('Name',['Sensitivity index of w over patch #' num2str(patch.number) ' for random variables #' num2str(alpha)])
        % set(gcf,'Name',['Sensitivity index of w over patch #' num2str(patch.number) ' for random variables #' num2str(alpha)])
    end
    clf
    v = variance(w{patch.number});
    s = varianceConditionalExpectation(w{patch.number},alpha)'./max(v);
    plot_sol(patch.S,s,varargin{:});
    colormap(p.Results.colormap)
    if p.Results.colorbar
        colorbar
    end
    set(gca,'FontSize',p.Results.FontSize)
end

end