function plot_sol(M,q,varargin)
% function plot_sol(M,q,varargin)

q = unfreevector(M,q);
options = patchoptions(getindim(M),varargin{:});
options = setcharin('ampl',options,getcharin('ampl',varargin,1));
if ischarin('sigma',varargin)
    options = [{'sigma' , getcharin('sigma',varargin)},options];
    options = setcharin('EdgeColor',options,getcharin('EdgeColor',varargin,'none'));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ischarin('epsilon',varargin)
    options = [{'epsilon' , getcharin('epsilon',varargin)},options];
    options = setcharin('EdgeColor',options,getcharin('EdgeColor',varargin,'none'));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
listegroup = getcharin('selgroup',varargin,1:getnbgroupelem(M));
for p=listegroup
    elem = getgroupelem(M,p);
    plot_lssol(elem,getnode(M),q,M.ls,options{:});
end
axis image
axis off

