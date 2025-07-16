function varargout = plot(u,M,varargin)
% function varargout = plot(u,M,'propertyname',propertyvalue,...)
% affichage d'un champ par element
% moyenne eventuellement si plusieurs points de gauss
% u : FEELEMFIELD
% M : MODEL
% propertyname       propertyvalue        effet
% 'selgroup'         liste                affiche les groupes d'elements de liste
% 'selelem'          liste                affiche les elements liste
% 'compo'            nom                  nom de la composante a afficher
% 'surface'          -----                pour un champ 2D , visualisation en 3D (marche pour un storage 'node')

if israndom(u)
    error('le FEELEMFIELD est aleatoire')
end

issurf = ischarin('surface',varargin);
varargin = delonlycharin('surface',varargin);

dim = getindim(M);
options = patchoptions(dim,'noedges',varargin{:});

numview = getcharin('view',varargin);
up_vector = getcharin('camup',varargin);
camera_position = getcharin('campos',varargin);
listegroup = getcharin('selgroup',varargin,1:M.nbgroupelem);
compo = getcharin('compo',varargin);
Helem = zeros(1,length(listegroup));
for p=listegroup
    elem = M.groupelem{p};
    if ischarin('selelem',varargin)
        listeelem =  getcharin('selelem',varargin,getnumelem(M)');
        elem = getelem(elem,listeelem,'global');
    end
    numlocal = getpos(M.groupelem{p},getnumber(elem));
    
    if ~strcmp(getlstype(elem),'out') || ischarin('out',varargin)
        if strcmp(u.type,'scalar')
            rep = 1;
            % ddlname = u.ddl;
        else
            ddlchoice = eval(['@get' u.type]);
            ddlchoice = ddlchoice(elem);
            if ~isempty(compo)
                [rep,repddlcompo] = findddl(ddlchoice,compo);
                if isempty(rep)
                    fprintf('Choisir une composantes parmi les suivantes :\n')
                    disp(get(ddlchoice,'ddl'))
                    error(' ')
                end
                % ddlname = getddlname(ddlchoice,rep);
            else
                rep = 1;
                % ddlname = getddlname(ddlchoice,rep) ;
            end
        end
        
        switch u.storage
            case {'center','gauss'}
                if strcmp(u.type,'scalar') && isa(u.value{p},'double')
                    value = reshape(u.value{p}(numlocal),[getnbelem(elem),1]);
                else
                    value = reshape(mean(u.value{p}(rep,:,numlocal,:),4),[getnbelem(elem),1]);
                end
                options = setcharin('FaceColor',options,'flat');
                
                if issurf && M.dim==2
                    Helem(p) = surface(elem,M.node,'FaceVertexCData',double(value),options{:});
                else
                    Helem(p) = plot(elem,M.node,'FaceVertexCData',double(value),options{:});
                end
            case {'node'}
                value = reshape(u.value{p}(rep,:,:),[M.nbnode,1]);
                options = setcharin('FaceColor',options,'interp');
                
                if issurf && M.dim==2
                    Helem(p) = surface(elem,M.node,'FaceVertexCData',double(value),options{:});
                    if isempty(numview)
                        view(3)
                    end
                else
                    Helem(p) = plot(elem,M.node,'FaceVertexCData',double(value),options{:});
                end
        end
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

% title(ddlname)
if ~issurf
    axis image
end
axis off

if nargout>=1
    varargout{1} = Helem;
end
