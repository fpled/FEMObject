function plotErrorStdAllSolutions(glob,patches,interfaces,U,w,lambda,U_ref,w_ref,lambda_ref,varargin)
% function plotErrorStdAllSolutions(glob,patches,interfaces,U,w,lambda,U_ref,w_ref,lambda_ref,varargin)
% Display the relative error in standard deviation of multiscale solution u=(U,w),
% global solution U, local solution w and Lagrange multiplier lambda
% with respect to reference multiscale solution u_ref=(U_ref,w_ref),
% global solution U_ref, local solution w_ref and Lagrange multiplier lambda_ref
% glob: Global
% patches: Patches
% interfaces: Interfaces
% U: FunctionalBasisArray of global solution U
% w: FunctionalBasisArray of local solution w
% lambda: FunctionalBasisArray of Lagrange multiplier lambda
% U_ref: FunctionalBasisArray of reference global solution U
% w_ref: FunctionalBasisArray of reference local solution w
% lambda_ref: FunctionalBasisArray of reference Lagrange multiplier lambda

p = ImprovedInputParser;
addParameter(p,'colorbar',true,@(x) islogical(x) || ischar(x));
addParameter(p,'colormap','default',@(x) isnumeric(x) || ischar(x));
addParameter(p,'FontSize',16,@isscalar);
parse(p,varargin{:})

varargin = delcharin({'colorbar','colormap','FontSize'},varargin);
n = numel(patches);

if ischarin('displ',varargin)
    i = getcharin('displ',varargin);
    figure('Name',['Relative error in standard deviation of u_' num2str(i) '=(U_' num2str(i) ',w_' num2str(i) '), U_' num2str(i) ', w_' num2str(i) ', lambda_' num2str(i) ' w.r.t. u_ref_' num2str(i) '=(U_ref_' num2str(i) ',w_ref_' num2str(i) '), U_ref_' num2str(i) ', w_ref_' num2str(i) ', lambda_ref_' num2str(i)])
    % set(gcf,'Name',['Relative error in standard deviation of u_' num2str(i) '=(U_' num2str(i) ',w_' num2str(i) '), U_' num2str(i) ', w_' num2str(i) ', lambda_' num2str(i) ' w.r.t. u_ref_' num2str(i) '=(U_ref_' num2str(i) ',w_ref_' num2str(i) '), U_ref_' num2str(i) ', w_ref_' num2str(i) ', lambda_ref_' num2str(i)])
else
    figure('Name','Relative error in standard deviation of u=(U,w), U, w, lambda w.r.t. u_ref=(U_ref,w_ref), U_ref, w_ref, lambda_ref')
    % set(gcf,'Name','Relative error in standard deviation of u=(U,w), U, w, lambda w.r.t. u_ref=(U_ref,w_ref), U_ref, w_ref, lambda_ref')
end
clf

if U.sz(1)==glob.S.nbddl
    P_out = calcProjection(glob,'free',false);
else
    P_out = glob.P_out;
end

subplot(1+n,2,1)
err_out = abs(std(U*P_out')-std(U_ref))'./max(std(U_ref));
err_out = unfreevector(glob.S_out,err_out)-calc_init_dirichlet(glob.S_out);
plot_sol(glob.S_out,err_out,varargin{:});
err_patch = cell(1,n);
for k=1:n
    patch = patches.patches{k};
    interface = interfaces.interfaces{k};
    err_patch{k} = abs(std(w{patch.number})-std(w_ref{patch.number}))'./max(std(w{patch.number}));
    plot_sol(patch.S,err_patch{k},varargin{:});
    if ~ischarin('sigma',varargin) && ~ischarin('epsilon',varargin) && ~ischarin('energyint',varargin)
        err_interface = abs(std(w{patch.number}*interface.P_patch')-std(w_ref{patch.number}*interface.P_patch'))'./max(std(w{patch.number}*interface.P_patch'));
        plot_sol(interface.S,err_interface,'FaceColor','none','EdgeColor','k',varargin{:});
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
plot_sol(glob.S_out,err_out,varargin{:});
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
    plot_sol(patch.S,err_patch{k},varargin{:});
    colormap(p.Results.colormap)
    if p.Results.colorbar
        colorbar
    end
    % axis(ax)
    % caxis(cax)
    set(gca,'FontSize',p.Results.FontSize)
    
    subplot(1+n,2,2*(k+1))
    interface = interfaces.interfaces{k};
    err_interface = abs(std(lambda{interface.number})-std(lambda_ref{interface.number}))'./max(std(lambda{interface.number}));
    plot_sol(interface.S,err_interface,'EdgeColor','interp',varargin{:});
    colormap(p.Results.colormap)
    if p.Results.colorbar
        colorbar
    end
    % axis(ax)
    % caxis(cax)
    set(gca,'FontSize',p.Results.FontSize)
end

end
