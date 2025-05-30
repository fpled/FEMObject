function plotClosedSobolIndicesAllSolutions(glob,patches,interfaces,U,w,lambda,alpha,varargin)
% function plotClosedSobolIndicesAllSolutions(glob,patches,interfaces,U,w,lambda,alpha,varargin)
% Display Closed Sobol indices of multiscale solution u=(U,w),
% global solution U, local solution w and Lagrange multiplier lambda
% associated with the group of variables alpha in {1,..,d} with d = min(ndims(U),ndims(w),ndims(lambda))
% glob: Global
% patches: Patches
% interfaces: Interfaces
% U: FunctionalBasisArray of global solution U
% w: FunctionalBasisArray of local solution w
% lambda: FunctionalBasisArray of Lagrange multiplier lambda
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
n = numel(patches);

if ischarin('displ',varargin)
    i = getcharin('displ',varargin);
    figure('Name',['Closed Sobol index of u_' num2str(i) '=(U_' num2str(i) ',w_' num2str(i) '), U_' num2str(i) ', w_' num2str(i) ', lambda_' num2str(i) ' for random variables #' num2str(alpha)])
    % set(gcf,'Name',['Closed Sobol index of u_' num2str(i) '=(U_' num2str(i) ',w_' num2str(i) '), U_' num2str(i) ', w_' num2str(i) ', lambda_' num2str(i) ' for random variables #' num2str(alpha)])
else
    figure('Name',['Closed Sobol index of u=(U,w), U, w, lambda for random variables #' num2str(alpha)])
    % set(gcf,'Name',['Closed Sobol index of u=(U,w), U, w, lambda for random variables #' num2str(alpha)])
end
clf

if U.sz(1)==glob.S.nbddl
    P_out = calcProjection(glob,'free',false);
else
    P_out = glob.P_out;
end

subplot(1+n,2,1)
dU = ndims(U);
sU_out = SensitivityAnalysis.closedSobolIndices(U*P_out',alpha,dU)';
sU_out = unfreevector(glob.S_out,sU_out)-calc_init_dirichlet(glob.S_out);
plot_sol(glob.S_out,sU_out,varargin{:});
for k=1:n
    patch = patches.patches{k};
    interface = interfaces.interfaces{k};
    dw = ndims(w{patch.number});
    plot_sol(patch.S,SensitivityAnalysis.closedSobolIndices(w{patch.number},alpha,dw)',varargin{:});
    if ~ischarin('sigma',varargin) && ~ischarin('epsilon',varargin) && ~ischarin('energyint',varargin)
        plot_sol(interface.S,SensitivityAnalysis.closedSobolIndices(w{patch.number}*interface.P_patch',alpha,dw)','FaceColor','none','EdgeColor','k',varargin{:});
    end
end
colormap(p.Results.colormap)
if p.Results.colorbar
    colorbar
end
ax=axis;
cax=caxis;
set(gca,'FontSize',p.Results.FontSize)

subplot(1+n,2,2)
sU = SensitivityAnalysis.closedSobolIndices(U,alpha,dU)';
sU = unfreevector(glob.S,sU)-calc_init_dirichlet(glob.S);
plot_sol(glob.S,sU,varargin{:});
colormap(p.Results.colormap)
if p.Results.colorbar
    colorbar
end
axis(ax)
% caxis(cax)
set(gca,'FontSize',p.Results.FontSize)

for k=1:n
    subplot(1+n,2,2*k+1)
    patch = patches.patches{k};
    dw = ndims(w{patch.number});
    plot_sol(patch.S,SensitivityAnalysis.closedSobolIndices(w{patch.number},alpha,dw)',varargin{:});
    colormap(p.Results.colormap)
    if p.Results.colorbar
        colorbar
    end
    % axis(ax)
    % caxis(cax)
    set(gca,'FontSize',p.Results.FontSize)
    
    subplot(1+n,2,2*(k+1))
    interface = interfaces.interfaces{k};
    dlambda = ndims(lambda{patch.number});
    plot_sol(interface.S,SensitivityAnalysis.closedSobolIndices(lambda{interface.number},alpha,dlambda)','EdgeColor','interp',varargin{:});
    colormap(p.Results.colormap)
    if p.Results.colorbar
        colorbar
    end
    % axis(ax)
    % caxis(cax)
    set(gca,'FontSize',p.Results.FontSize)
end

end
