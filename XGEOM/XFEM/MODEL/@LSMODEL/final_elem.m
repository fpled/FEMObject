function [M,I_elemtot,H]= final_elem(M,varargin)

% M= final_elem(M,PCls,option)
% finalise le MODEL M en créant les partitions stochastiques de chaque EF
% -> renumerotation des noeuds de 1 a M.nbnode
% -> calcul des connectivites des noeuds et elements
% -> creation des ddl

M = unique(M);
M = changenodenumber(M,[1:M.nbnode]);
M = changeelemnumber(M);
M = calc_connec(M);

% ELIMINATION DES NOEUDS SANS ELEMENT
repnodeelim = find(sum(M.connec.node2elem,1)==0);
if any(repnodeelim)

    M = removenode(M,getnumber(M.node,repnodeelim));
    M = calc_connec(M);
%%%%
 %M = changenodenumber(M,[1:M.nbnode]);
%%%%
end

% RENUMEROTATION POUR REDUIRE LA LARGEUR DE BANDE
if ischarin('renum',varargin)
 disp('Warning : renumerotation avec symrcm')
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

%GESTION DES LEVELSETS

M.ls=LEVELSETS(M.ls);
if length(M.ls)>0
[M,I_elemtot,H] = lssplitelem_elem(M,varargin{:});
varargout{1}=H;
%M = PCMODEL(M,'manudim',getnumber(H),H,varargin{:});
end

M = actualisematerials(M);


% CREATION DES DDL
M=createddlnode(M,varargin{:});

M = lsbloqueoutddl(M);

