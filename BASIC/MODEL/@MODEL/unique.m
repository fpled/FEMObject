function [M,nodeelim,nodereplace]=unique(M,fununique,display)
% [M,nodeelim,nodereplace]=unique(M,fununique,display)
% Pour trouver les noeuds doubles on regarde la distance des noeuds par rapport a deux
% noeuds distincts. Si les distances de plusieurs points aux 2 noeuds de
% reference sont egales, on ne conserve qu'un des noeuds sur cet ensemble.
if nargin<2
    display = false;
    if nargin==1
        fununique=@(x)unique(x);
    end
end

[~,repelim,repreplace]=fununique(POINT(M.node));

if any(repelim)
    nodeelim=getnumber(M.node,repelim);
    if display
        disp(['!!! elimination des double noeuds : ',...
            num2str(length(nodeelim)) ' noeuds supprimes'])
    end
    nodereplace = getnumber(M.node,repreplace);
    M = removenode(M,nodeelim,nodereplace);
    
    for p=1:M.nbgroupelem
        M.groupelem{p}=unique(M.groupelem{p});
        for j=1:p-1
            if strcmp(class(M.groupelem{j}),class(M.groupelem{p}))
                M.groupelem{p}=setdiff(M.groupelem{p},M.groupelem{j});
            end
        end
    end
end