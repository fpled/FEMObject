function varargout = plot(M,varargin)
% function [Helem] = plot(M,'propertyname',propertyvalue,...)
% M : MODEL
% propertyname       propertyvalue        effet
% 'selgroup'         liste                affiche les groupes d'elements de liste
% 'selelem'          liste                affiche les elements liste
% 'node'             ---                  affiche les noeuds avec des petits points
% 'numnode'          ---                  affichage des numeros des noeuds
% 'numelem'          ---                  affichage des numeros des elements
% 'mat'              ---                  affichage des materiaux
% 'group'            ---                  affichage des groupes d'elements
% 'color'            ---                  definit la couleur du maillage
% -> voir patchoptions pour la suite des parametres associes a la commande patch
%
% [Helem,Hnode] = function plot(M,'propertyname',propertyvalue,...)
% Helem : handle vers les patchs d'elements

if ~isa(M,'MODEL')
    if nargin>=2 && isa(varargin{1},'MODEL')
        % try
        if numel(M)==getnbelem(varargin{1})
            Helem = plot(FEELEMFIELD(M,varargin{:}),varargin{:});
        else
            Helem = plot(FENODEFIELD(unfreevector(varargin{1},M)),varargin{:});
        end
        
        % catch
        %     error('plot n''a pas marche')
        % end
    else
        error('plot n''est pas possible')
    end
else
    issurf = ischarin('surface',varargin);
    varargin = delonlycharin('surface',varargin);
    
    nodetext = ischarin('numnode',varargin);
    elemtext = ischarin('numelem',varargin);
    matplot = ischarin('mat',varargin);
    lsenrichplot = ischarin('lsenrich',varargin);
    groupplot = ischarin('group',varargin);
    listegroup = getcharin('selgroup',varargin,1:M.nbgroupelem);
    listeelem = getcharin('selelem',varargin,getnumelem(M));
    
    dim = getindim(M);
    options = patchoptions(dim,varargin{:});
    
    color = getcharin('color',varargin);
    if ~isempty(color) && isa(color,'char')
        options = setcharin('edgecolor',options,color);
    elseif ~isempty(color) && isa(color,'double')
        options = setcharin('facecolor',options,'flat');
        options = setcharin('facevertexcdata',options,color);
    end
    node = M.node;
    nodecoord = double(getcoord(node));
    
    Helem = [];
    Hnode = [];
    for p=listegroup
        if ischarin('selelem',varargin)
            if ischarin('local',varargin)
                elem = getelem(M.groupelem{p},listeelem,'local');
            else
                elem = getelem(M.groupelem{p},listeelem,'global');
            end
        else
            elem =  M.groupelem{p};
        end
        optionselem = options;
        if getdim(elem)==1
            optionselem = delcharin('facelighting',optionselem);
            optionselem = delcharin('edgelighting',optionselem);
        end
        if getnbelem(elem)>0 && ~(strcmp(getlstype(elem),'out') && ischarin('noout',varargin))
            
            if matplot || groupplot || lsenrichplot
                if matplot
                    colelem = getmaterialnumber(elem);
                elseif groupplot
                    colelem = p;
                elseif lsenrichplot
                    colelem = getlsenrich(elem);
                end
                if getdim(elem)==1
                    optionselem = setcharin('edgecolor',optionselem,getfacecolor(colelem));
                else
                    optionselem = setcharin('facecolor',optionselem,'flat');
                    optionselem = setcharin('facevertexcdata',optionselem,colelem);
                end
            end
            
            Helemtemp = plot(elem,node,optionselem{:});
            Helem = [Helem,Helemtemp];
            
            if elemtext==1
                plotnumber(elem,node,'color','r',varargin{:});
            end
            
            if  getdim(elem)>1 && (matplot || groupplot || lsenrichplot)
                xplot = calc_midpoint(elem,node);
                plottext(xplot,num2str(colelem),'color','w');
            end
            
        end
    end
    
    if nodetext==1
        plotnumber(node,varargin{:});
    end
    
    axis image
    axis off
    
    numview = getcharin('view',varargin);
    up_vector = getcharin('camup',varargin);
    camera_position = getcharin('campos',varargin);
    if ~isempty(numview)
        view(numview)
    elseif dim==3 || issurf
        view(3)
    end
    if ~isempty(up_vector)
        camup(up_vector)
    end
    if ~isempty(camera_position)
        campos(camera_position);
    end
end

if nargout>=1
    varargout{1} = Helem;
end
