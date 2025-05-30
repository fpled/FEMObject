function [numelem,numnode,numgroupelem] = findeleminnode(M,node)
% function [numelem,numnode,numgroupelem] = findeleminnode(M,node)
% renvoie l'ensemble numelem des elements de M dont tous les noeuds appartiennent a l'ensemble node,
% l'ensemble numnode des noeuds de M connectes a ces elements
% et l'ensemble numgroupelem des groupes d'elements correspondants

if isa(node,'NODE')
    node = getnumber(node);
end
numelem = [];
numnode = [];
numgroupelem = [];
for j=1:M.nbgroupelem
    [numelemtemp,numntemp] = findeleminnode(M.groupelem{j},node);
    numelem = [numelem;numelemtemp];
    numnode = union(numnode,numntemp);
    if ~isempty(numelemtemp)
        numgroupelem = [numgroupelem,j];
    end
end
