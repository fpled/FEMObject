function subgauss = elem_subgauss(elem,p,m,choix,xnodereel)
%function subgauss = elem_subgauss(elem,p,m,choix,xnode)
% p : ordre
% m : nombre de sous-decoupages
% si choix=2 : calcul des poids et coordonnees des sous-points de gauss (par defaut)
% si choix=1 : calcul des poids
%
% xnodereel (optionnel) : coordonnees des noeuds de l'elements
%     dans ce cas points de gauss dans le repere reel

if nargin<=3 || isempty(choix)
    choix = 2;
end

if nargin<=4
    xnodereel = nodelocalcoordtet4();
end

subgauss.w = [];
subgauss.coord = zeros(0,3);

subgauss = recursivesubgauss(subgauss,p,m,0,xnodereel,choix);

subgauss.nbgauss = numel(subgauss.w);

function subgauss = recursivesubgauss(subgauss,p,m,i,xnode,choix)

if i<m
    [subconnec,xcenter] = splitpolyhedron(1:4,xnode,5:10);
    
    xnodeplus = [xnode;xcenter];
    for k=1:size(subconnec,1)
        xnode=xnodeplus(subconnec(k,:),:);
        subgauss = recursivesubgauss(subgauss,p,m,i+1,xnode,choix);
    end
    
    %col = {'y','m','g','r','b','k','g','c'};
    %for k=1:size(subconnec,1)
    %    plot(TET4(NODE(xnodeplus),1,subconnec(k,:)),NODE(xnodeplus),'FaceColor',col{k},'FaceLighting','gouraud','FaceAlpha',.5);
    %end
    
else
    tet4gauss = gausstet4(p);
    w = tet4gauss.w * detJtet4(xnode);
    subgauss.w = [subgauss.w;w(:)] ;
    if choix==2
        c = Ntet4(tet4gauss.coord)*xnode;
        subgauss.coord = [subgauss.coord;c] ;
    end
end

return

function [subconnec,xplus] = splitpolyhedron(connec,xnode,numplus)
numplus = 5:10;
xplus = [1/2,0,0;1/2,1/2,0;0,1/2,0;1/2,0,1/2;0,0,1/2;0,1/2,1/2];
xplus = Ntet4(xplus)*xnode;
subconnec = [1,5,7,9;...
    5,2,6,8;...
    7,6,3,10;...
    9,8,10,4;...
    7,8,6,10;...
    5,6,7,8;...
    7,8,10,9;...
    5,7,9,8];
return

function [subconnec,xplus] = splitpolyhedron1(connec,xnode,numplus)
numplus = 5:7;
xplus = [1/2,1/2,0;0,1/2,1/2;1/2,0,1/2];
xplus = Ntet4(xplus)*xnode;
subconnec = [connec(1),connec(2),numplus(1),numplus(3);...
    connec(1),numplus(1),connec(3),numplus(2);...
    connec(1),numplus(3),numplus(2),connec(4);...
    connec(1),numplus(1),numplus(2),numplus(3)];
return
