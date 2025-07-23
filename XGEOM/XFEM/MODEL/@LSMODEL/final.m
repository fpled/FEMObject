function [M,varargout] = final(M,varargin)
% M = final(M,option)
% finalise le LSMODEL M
% -> renumerotation des noeuds de 1 a M.nbnode
% -> calcul des connectivites des noeuds et elements
% -> creation des ddl

M = unique(M);
% M = changenodenumber(M,[1:getnbnode(M)]);
M = sortnodenumber(M);
M = changeelemnumber(M);
M = calc_connec(M);

% ELIMINATION DES NOEUDS SANS ELEMENT
connec = getconnec(M);
repnodeelim = find(sum(connec.node2elem,1)==0);
if any(repnodeelim)
    M = removenode(M,getnumber(getnode(M),repnodeelim));
    M = calc_connec(M);
end

% RENUMEROTATION POUR REDUIRE LA LARGEUR DE BANDE
if ischarin('renum',varargin)
    disp('Warning : renumerotation avec symrcm')
    % essayer le cas 2D domaine [0,1]^2, PLANE CONT avec poisson 0.3
    nodenum = getnumber(getnode(M));
    connec = getconnec(M);
    num = symrcm(connec.node2node); % reduire la largeur de bande
    [r,renum]=sort(num);
    M = changenodenumber(M,nodenum,renum);
    connec.node2node = connec.node2node(num,num);
    connec.node2elem = connec.node2elem(:,num);
    connec.elem2node = connec.elem2node(num,:);
end
varargin = delonlycharin({'norenum','renum'},varargin);

%GESTION DES LEVELSETS

if getnblevelsets(M)>0
    M = lssplitelem(M,varargin{:});
    
    %M = applyfunctiontofaces(M,@lssplitelem,varargin{:});
    
    if ~ischarin('norenumelem',varargin)
        M = changeelemnumber(M);
        M = calc_connec(M);
    end
end

try
    M = actualisematerials(M);
end

M = lsenrich(M,varargin{:});

M = createddlnode(M,varargin{:});
M = changeelemnumber(M);
M = calc_connec(M);

if ~ischarin('nobloqueout',varargin)
    M = lsbloqueoutddl(M);
end

