function plotClosedSobolIndicesLocalSolution(patches,w,alpha,varargin)
% function plotClosedSobolIndicesLocalSolution(patches,w,alpha,varargin)
% Display the Closed Sobol indices of local solution w associated with the
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
        figure('Name',['Closed Sobol index of w_' num2str(i) ' over patches #' num2str([numbers{:}]) ' for random variables #' num2str(alpha)])
        % set(gcf,'Name',['Closed Sobol index of w_' num2str(i) ' over patches #' num2str([numbers{:}]) ' for random variables #' num2str(alpha)])
    else
        figure('Name',['Closed Sobol index of w over patches #' num2str([numbers{:}]) ' for random variables #' num2str(alpha)])
        % set(gcf,'Name',['Closed Sobol index of w over patches #' num2str([numbers{:}] ' for random variables #' num2str(alpha)])
    end
    clf
    n = numel(patches);
    for k=1:n
        patch = patches.patches{k};
        d = ndims(w{patch.number});
        plot_sol(patch.S,SensitivityAnalysis.closedSobolIndices(w{patch.number},alpha,d)',varargin{:});
        colormap(p.Results.colormap)
        if p.Results.colorbar
            colorbar
        end
        set(gca,'FontSize',p.Results.FontSize)
    end
elseif isa(patches,'Patch')
    patch = patches;
    d = ndims(w{patch.number});
    if ischarin('displ',varargin)
        i = getcharin('displ',varargin);
        figure('Name',['Closed Sobol index of w_' num2str(i) ' over patch #' num2str(patch.number) ' for random variables #' num2str(alpha)])
        % set(gcf,'Name',['Closed Sobol index of w_' num2str(i) ' over patch #' num2str(patch.number) ' for random variables #' num2str(alpha)])
    else
        figure('Name',['Closed Sobol index of w over patch #' num2str(patch.number) ' for random variables #' num2str(alpha)])
        % set(gcf,'Name',['Closed Sobol index of w over patch #' num2str(patch.number) ' for random variables #' num2str(alpha)])
    end
    clf
    plot_sol(patch.S,SensitivityAnalysis.closedSobolIndices(w{patch.number},alpha,d)',varargin{:});
    colormap(p.Results.colormap)
    if p.Results.colorbar
        colorbar
    end
    set(gca,'FontSize',p.Results.FontSize)
end

end
