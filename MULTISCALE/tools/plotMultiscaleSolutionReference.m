function plotMultiscaleSolutionReference(glob,patches,interfaces,U_ref,w_ref,varargin)
% function plotMultiscaleSolutionReference(glob,patches,interfaces,U_ref,w_ref,varargin)
% Display reference multiscale solution u_ref=(U_ref,w_ref)
% glob: Global or GlobalOutside
% patches: Patches
% interfaces: Interfaces
% U_ref: reference global solution U
% w_ref: reference local solution w

p = ImprovedInputParser;
addParameter(p,'colorbar',true,@(x) islogical(x) || ischar(x));
addParameter(p,'colormap','default',@(x) isnumeric(x) || ischar(x));
addParameter(p,'FontSize',16,@isscalar);
parse(p,varargin{:})

varargin = delcharin({'colorbar','colormap','FontSize'},varargin);
n = numel(patches);

if ischarin('sigma',varargin)
    i = getcharin('sigma',varargin);
    figure('Name',['Reference multiscale solution sig_u_ref_' num2str(i) '=(sig_U_ref_' num2str(i) ',sig_w_ref_' num2str(i) ') over complementary subdomain and patches'])
    % set(gcf,'Name',['Reference multiscale solution sig_u_ref_' num2str(i) '=(sig_U_ref_' num2str(i) ',sig_w_ref_' num2str(i) ') over complementary subdomain and patches'])
elseif ischarin('epsilon',varargin)
    i = getcharin('epsilon',varargin);
    figure('Name',['Reference multiscale solution eps_u_ref_' num2str(i) '=(eps_U_ref_' num2str(i) ',eps_w_ref_' num2str(i) ') over complementary subdomain and patches'])
    % set(gcf,'Name',['Reference multiscale solution eps_u_ref_' num2str(i) '=(eps_U_ref_' num2str(i) ',eps_w_ref_' num2str(i) ') over complementary subdomain and patches'])
elseif ischarin('energyint',varargin)
    figure('Name',['Reference multiscale solution H_u_ref=(H_U_ref,H_w_ref) over complementary subdomain and patches'])
    % set(gcf,'Name',['Reference multiscale solution H_u_ref=(H_U_ref,H_w_ref) over complementary subdomain and patches'])
elseif ischarin('displ',varargin)
    i = getcharin('displ',varargin);
    figure('Name',['Reference multiscale solution u_ref_' num2str(i) '=(U_ref_' num2str(i) ',w_ref_' num2str(i) ') over complementary subdomain and patches'])
    % set(gcf,'Name',['Reference multiscale solution u_ref_' num2str(i) '=(U_ref_' num2str(i) ',w_ref_' num2str(i) ') over complementary subdomain and patches'])
elseif ischarin('rotation',varargin)
    i = getcharin('rotation',varargin);
    figure('Name',['Reference multiscale solution rot_u_ref_' num2str(i) '=(rot_U_ref_' num2str(i) ',rot_w_ref_' num2str(i) ') over complementary subdomain and patches'])
    % set(gcf,'Name',['Reference multiscale solution rot_u_ref_' num2str(i) '=(rot_U_ref_' num2str(i) ',rot_w_ref_' num2str(i) ') over complementary subdomain and patches'])
else
    figure('Name','Reference multiscale solution u_ref=(U_ref,w_ref) over complementary subdomain and patches')
    % set(gcf,'Name','Reference multiscale solution u_ref=(U_ref,w_ref) over complementary subdomain and patches')
end
clf

if isa(glob,'Global')
    S_out = glob.S_out;
elseif isa(glob,'GlobalOutside')
    S_out = glob.S;
end

plot_sol(S_out,U_ref,varargin{:});
for k=1:n
    patch = patches.patches{k};
    interface = interfaces.interfaces{k};
    plot_sol(patch.S,w_ref{patch.number},varargin{:});
    if ~ischarin('sigma',varargin) && ~ischarin('epsilon',varargin) && ~ischarin('energyint',varargin)
        plot_sol(interface.S,interface.P_patch*w_ref{patch.number},'FaceColor','none','EdgeColor','k',varargin{:});
    end
end
colormap(p.Results.colormap)
if p.Results.colorbar
    colorbar
end
set(gca,'FontSize',p.Results.FontSize)

end
