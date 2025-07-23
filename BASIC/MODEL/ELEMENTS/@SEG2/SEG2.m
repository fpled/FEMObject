function elem = SEG2(node,numelem,connec,varargin)
% function elem = SEG2(node,numelem,connec,'material',mat,'option',option,syscoord,'parent',parent)
% node : objet de type NODE contenant les noeuds de l'element
% connec : table de connectivite
% numelem : numero des elements
% option : 'BORD' si c'est un element de bord

if nargin==0    
    elemp = ELEMENTGEOM(1);
    eleme = ELEMENT();
    
    elem = struct();
    elem = class(elem,'SEG2',elemp,eleme);
else    
    if nargin<2
        numelem = 1;
    end

    if nargin<3 && numelem==1
        connec = [1,2];
    elseif nargin<3
        error('rentrer la connectivite')
    end

    if isa(node,'NODE')
        connec = reshape(connec,length(numelem),2);
        P1 = POINT(getnode(node,connec(:,1)));
        P2 = POINT(getnode(node,connec(:,2)));
    elseif isa(node,'POINT')
        P2 = node(1,:);
        P1 = node(2,:);
    elseif isa(node,'double') || isa(node,'MYDOUBLEND')
        P2 = POINT(node(1,:));
        P1 = POINT(node(2,:));
    end
    eX = P2-P1;
    param.L = distance(P1,P2);

    eX = normalize(eX);
    syscoordlocal = CARTESIAN1D(eX);
    if isempty(connec)
        syscoord = CARTESIAN3D();
    else
        switch getindim(eX)
        case 1
            syscoord = CARTESIAN1D(eX);        
        case 2
            eY = rot2D(eX,pi/2);
            syscoord = CARTESIAN2D(eX,eY); 
        case 3
            if ischarin('param',varargin) && isa(getcharin('param',varargin),'VECTOR')
                eY = getcharin('param',varargin);
                eZ = cross(eX,eY);
                eZ = normalize(eZ);
                eY = cross(eZ,eX);
            else
                [eY,eZ] = planortho(eX);
            end
            syscoord = CARTESIAN3D(eX,eY,eZ); 
        end
    end
    
    elemp = ELEMENTGEOM(1,node,numelem,connec,syscoordlocal,syscoord);
    eleme = ELEMENT(varargin{:});
    eleme = setparam(eleme,'L',param.L);
    
    elem = struct();
    elem = class(elem,'SEG2',elemp,eleme);
end
