function [M,varargout] = final(M,varargin)
% M = final(M,option)
% finalise le MODEL M
% -> renumerotation des noeuds de 1 a M.nbnode
% -> calcul des connectivites des noeuds et elements
% -> creation des ddl

if ~ischarin('duplicate',varargin)
    M = unique(M);
end
varargin = delonlycharin('duplicate',varargin);
% M = changenodenumber(M,int32([1:M.nbnode]));
M = sortnodenumber(M);
M = changeelemnumber(M);
M = calc_connec(M);

% ELIMINATION DES NOEUDS SANS ELEMENT
repnodeelim = find(sum(M.connec.node2elem,1)==0);
if any(repnodeelim)
    M = removenode(M,getnumber(M.node,repnodeelim));
    M = calc_connec(M);
end

% RENUMEROTATION POUR REDUIRE LA LARGEUR DE BANDE
if ischarin('renum',varargin)
    warning('Renumerotation avec symrcm');
    % essayer le cas 2D domaine [0,1]^2, PLANE CONT avec poisson 0.3
    nodenum = getnumber(M.node);
    num = symrcm(M.connec.node2node); % reduire la largeur de bande
    [r,renum]=sort(num);
    M = changenodenumber(M,nodenum,renum);
    M.connec.node2node = M.connec.node2node(num,num);
    M.connec.node2elem = M.connec.node2elem(:,num);
    M.connec.elem2node = M.connec.elem2node(num,:);
end
varargin = delonlycharin({'norenum','renum'},varargin);

try
    M = actualisematerials(M);
end
% CREATION DES DDL

M = createddlnode(M,varargin{:});

