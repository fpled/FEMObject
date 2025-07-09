function contour(u,M,varargin)
% function contour(u,M,varargin)
% affichage d'un champ defini aux noeuds
% u : FENODEFIELD
% M : MODEL
% argi = 'edges' affichage des bords des elements en 2D
% argi = 'surface' pour un champ 2D , visualisation en 3D
% argi = 'surfacemesh' pour un champ 2D , visualisation en 3D et affichage
% du maillage dans le plan z=0
% argi = 'nofaces' affichage des bords des elements en 2D

facecolor = getcharin('facecolor',varargin);
if isempty(facecolor)
    if ischarin('nofaces',varargin)
        facecolor = 'none';
    else
        facecolor = 'interp';
    end
end

edgecolor = getcharin('edgecolor',varargin);
if isempty(edgecolor)
    if ischarin('edges',varargin)
        edgecolor = 'k';
    else
        edgecolor = 'none';
    end
end

numview = getcharin('view',varargin);
up_vector = getcharin('camup',varargin);
camera_position = getcharin('campos',varargin);

for p=1:M.nbgroupelem
    if (ischarin('surface',varargin) || ischarin('surfacemesh',varargin)) && M.dim==2
        surface(M.groupelem{p},M.node,'facecolor',facecolor,'edgecolor',edgecolor,'facevertexcdata',u.value);
        
        if ischarin('surfacemesh',varargin)
            plot(M.groupelem{p},M.node,'facecolor','none','edgecolor','k','facevertexcdata',u.value);
        end
        axis image
        axis square
        if isempty(numview)
            view(3)
        end
    else
        plot(M.groupelem{p},M.node,'facecolor',facecolor,'edgecolor',edgecolor,'facevertexcdata',u.value);
        axis image
    end
end

axis off

if ~isempty(numview)
    view(numview)
elseif M.indim==3
    view(3)
end
if ~isempty(up_vector)
    camup(up_vector)
end
if ~isempty(camera_position)
    campos(camera_position);
end

