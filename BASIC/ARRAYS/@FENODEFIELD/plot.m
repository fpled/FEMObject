function varargout = plot(u,M,varargin)
% function varargout = plot(u,M,'propertyname',propertyvalue,...)
% affichage d'un champ defini aux noeuds
% u : FENODEFIELD
% M : MODEL
% propertyname       propertyvalue        effet
% 'selgroup'         liste                affiche les groupes d'elements de liste
% 'selelem'          liste                affiche les elements liste
% 'surface'          -----                pour un champ 2D , visualisation en 3D
% 'surfacemesh'      -----                comme 'surface' plus affichage du maillage

if israndom(u)
    error('le FENODEFIELD est aleatoire')
end
u = double(u);
if length(u)~=getnbnode(M)
    u = unfreevector(M,u);
end

dim = getindim(M);
options = patchoptions(dim,'noedges',varargin{:});
% options = setcharin('EdgeColor',options,'none');
if strcmp(getcharin('FaceColor',options),'none')
    options = setcharin('FaceColor',options,'interp');
end

numview = getcharin('view',varargin);
up_vector = getcharin('camup',varargin);
camera_position = getcharin('campos',varargin);
listegroup = getcharin('selgroup',varargin,1:M.nbgroupelem);
Helem = zeros(1,length(listegroup));
for p=listegroup
    elem = M.groupelem{p};
    if ischarin('selelem',varargin)
        listeelem =  getcharin('selelem',varargin,getnumelem(M)');
        if ischarin('local',varargin)
            elem = getelem(elem,listeelem,'local');
        else
            elem = getelem(elem,listeelem,'global');
        end
    end
    
    if dim==1 && ischarin('courbe',varargin)
        x = double(getcoord(M.node));
        options = delonlycharin('courbe',varargin);
        Helem(p) = plot(x,u,options{:});
        
    elseif (ischarin('surface',varargin) || ischarin('surfacemesh',varargin)) && dim==2
        options = delonlycharin('surface',options);
        Helem(p) = surface(elem,M.node,'FaceVertexCData',u,options{:});
        
        if ischarin('surfacemesh',varargin)
            options = setcharin('EdgeColor',options,'k');
            options = setcharin('FaceColor',options,'none');
            options = delonlycharin('surfacemesh',options);
            plot(elem,M.node,'FaceVertexCData',u,options{:});
        end
        axis image
        if isempty(numview)
            view(3)
        end
    else
        Helem(p) = plot(elem,M.node,'FaceVertexCData',full(u),options{:});
        axis image
    end
end

if ~isempty(numview)
    view(numview)
elseif dim==3
    view(3)
end
if ~isempty(up_vector)
    camup(up_vector)
end
if ~isempty(camera_position)
    campos(camera_position);
end

if ~(dim==1 && ischarin('courbe',varargin))
    axis off
end

if nargout>=1
    varargout{1} = Helem;
end
