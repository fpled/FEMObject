function plotErrorMeanMultiscaleSolution(glob,patches,interfaces,U,w,U_ref,w_ref,varargin)
% function plotErrorMeanMultiscaleSolution(glob,patches,interfaces,U,w,U_ref,w_ref,varargin)
% Display the relative error in the mean (mathematical expectation) of multiscale solution u=(U,w)
% with respect to reference multiscale solution u_ref=(U_ref,w_ref)
% glob: Global
% patches: Patches
% U: FunctionalBasisArray of global solution U
% w: FunctionalBasisArray of local solution w
% U_ref: FunctionalBasisArray of reference global solution U
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
    figure('Name',['Relative error in mean of u_' num2str(i) '=(U_' num2str(i) ',w_' num2str(i) ') w.r.t. u_ref_' num2str(i) '=(U_ref_' num2str(i) ',w_ref_' num2str(i) ')'])
    % set(gcf,'Name',['Relative error in mean of u_' num2str(i) '=(U_' num2str(i) ',w_' num2str(i) ') w.r.t. u_ref_' num2str(i) '=(U_ref_' num2str(i) ',w_ref_' num2str(i) ')'])
else
    figure('Name','Relative error in mean of u=(U,w) w.r.t. u_ref=(U_ref,w_ref)')
    % set(gcf,'Name','Relative error in mean of u=(U,w) w.r.t. u_ref=(U_ref,w_ref)')
end
clf

if U.sz(1)==glob.S.nbddl
    P_out = calcProjection(glob,'free',false);
else
    P_out = glob.P_out;
end

err_out = abs(mean(U*P_out')-mean(U_ref))'./max(mean(U_ref));
err_out = unfreevector(glob.S_out,err_out)-calc_init_dirichlet(glob.S_out);
plot_sol(glob.S_out,err_out,varargin{:});
for k=1:n
    patch = patches.patches{k};
    interface = interfaces.interfaces{k};
    err_patch = abs(mean(w{patch.number})-mean(w_ref{patch.number}))'./max(mean(w{patch.number}));
    plot_sol(patch.S,err_patch,varargin{:});
    if ~ischarin('sigma',varargin) && ~ischarin('epsilon',varargin) && ~ischarin('energyint',varargin)
        err_interface = abs(mean(w{patch.number}*interface.P_patch')-mean(w_ref{patch.number}*interface.P_patch'))'./max(mean(w{patch.number}*interface.P_patch'));
        plot_sol(interface.S,err_interface,'FaceColor','none','EdgeColor','k',varargin{:});
    end
end
colormap(p.Results.colormap)
if p.Results.colorbar
    colorbar
end
set(gca,'FontSize',p.Results.FontSize)

end
