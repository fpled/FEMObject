function plot_enrichfun(M,lsnum,varargin)
% function plot_enrichfun(M,lsnum,varargin)
% lsnum : numero de la levelset concernee

issurf = ischarin('surface',varargin);
iscourbe = ischarin('courbe',varargin);
varargin = delonlycharin('surface',varargin);
dim = getindim(M);
options = patchoptions(dim,varargin{:});
options = setcharin('FaceColor',options,'interp');
options = setcharin('EdgeColor',options,'none');

if issurf
    options = [options , {'surface'}];
end
if iscourbe
    options = [options , {'courbe'}];
end
listegroup =  getcharin('selgroup',varargin,1:M.nbgroupelem);
for p=listegroup
    elem = getgroupelem(M,p);
    plot_enrichfun(elem,getnode(M),M.ls,lsnum,options{:});
end

if ~iscourbe
    axis image
    axis off
    if issurf
        axis square
        view(3)
    end
end
