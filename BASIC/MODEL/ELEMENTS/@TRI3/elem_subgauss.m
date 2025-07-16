function subgauss = elem_subgauss(elem,p,m,choix,xnodereel)
%function subgauss = elem_subgauss(elem,p,m,choix,xnodereel)
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
    xnodereel = nodelocalcoordtri3();
end

subgauss.w = [];
subgauss.coord = zeros(0,2);

subgauss = recursivesubgauss(subgauss,p,m,0,xnodereel,choix);

subgauss.nbgauss = numel(subgauss.w);

function subgauss = recursivesubgauss(subgauss,p,m,i,xnode,choix)

if i<m
    [subconnec,xcenter]=splitpoly2(1:3,xnode,4);
    xnodeplus = [xnode;xcenter];
    for k=1:size(subconnec,1)
        xnode = xnodeplus(subconnec(k,:),:);
        subgauss = recursivesubgauss(subgauss,p,m,i+1,xnode,choix);
    end
else
    tri3gauss = gausstri3(p);
    w = tri3gauss.w * detJtri3(xnode);
    subgauss.w = [subgauss.w;w(:)] ;
    if choix==2
        c = Ntri3(tri3gauss.coord)*xnode;
        subgauss.coord = [subgauss.coord;c] ;
    end
    % plot(TRI3(NODE(xnode),1,1:3),NODE(xnode),'FaceColor','none');
end



return


function [subconnec,xplus] = splitpoly2(connec,xnode,numplus)
numplus = 4:6;
xplus = [1/2,0;1/2,1/2;0,1/2];
xplus = Ntri3(xplus)*xnode;
subconnec = [connec(1),numplus(1),numplus(3);...
    connec(2),numplus(2),numplus(1);...
    connec(3),numplus(3),numplus(2);...
    numplus(1),numplus(2),numplus(3)];
return
