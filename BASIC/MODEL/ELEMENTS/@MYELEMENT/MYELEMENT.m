function elem = MYELEMENT(node,numelem,connec,mat,option,varargin)
% function elem = MYELEMENT(node,numelem,connec,mat,option,varargin)
%
% node : objet de type NODE contenant les noeuds de l'element
% matnum : numero du materiau
% option : ''
% connec : table de connectivite
% numelem : numero des elements
% eY (uniquement en 3D) : eX etant le vecteur tangent a l'element, 
%       eY est un vecteur permettant de definir le plan principal (eX,eY) de la poutre


elemg = SEG2(node,numelem,connec,mat,option);

eX = VECTOR(getsyscoordlocal(elemg));
param.L = diameter(elemg,node);

switch getindim(eX)
case 2
    eY = rot2D(eX,pi/2);
    param.syscoordlocalsect = CARTESIAN2D(eX,eY); 
case 3
    try
        eY = varargin{1};
        eZ = cross(eX,eY);
        eZ = normalize(eZ);
        eY = cross(eZ,eX);
    catch
        fprintf('WARNING : systeme de coordonnees local au segment defini arbitrairement\n')
        [eY,eZ] = planortho(eX);
    end
    param.syscoordlocalsect = CARTESIAN3D(eX,eY,eZ); 
end

elem.param = param ;

elem = class(elem,'MYELEMENT',elemg);
