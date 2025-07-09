function varargout = plot_sol(M,q,varargin)
% function varargout = plot_sol(M,q,varargin)

if ~isa(M,'MODEL')
    varargout = cell(1,nargout);
    [varargout{:}] = plot_sol(q,M,varargin{:});
    return
end

if israndom(q)
    error('la solution est aleatoire')
end
q = unfreevector(M,q);

dim = getdim(M);
options = patchoptions(dim,'noedges',varargin{:});
options = setcharin('ampl',options,getcharin('ampl',varargin,0));
if ischarin('sigma',varargin)
    options = [{'sigma',getcharin('sigma',varargin)},options];
elseif ischarin('epsilon',varargin)
    options = [{'epsilon',getcharin('epsilon',varargin)},options];
elseif ischarin('energyint',varargin)
    options = [{'energyint',getcharin('energyint',varargin)},options];
elseif ischarin('displ',varargin)
    options = [{'displ',getcharin('displ',varargin)} ,options];
elseif ischarin('rotation',varargin)
    options = [{'rotation',getcharin('rotation',varargin)} ,options];
end

listegroup = getcharin('selgroup',varargin,1:M.nbgroupelem);
Helem = zeros(1,length(listegroup));
for p=listegroup
    elem = M.groupelem{p};
    if ischarin('selelem',varargin)
        listeelem = getcharin('selelem',varargin,getnumelem(M)');
        if ischarin('local',varargin)
            elem = getelem(elem,listeelem,'local');
        else
            elem = getelem(elem,listeelem,'global');
        end
    end
    
    if dim==1 && ischarin('courbe',varargin)
        x = getcoord(M.node);
        options = delonlycharin('courbe',varargin);
        Helem(p) = plot(x,q,options{:});
    else
        Helem(p) = plot_sol(elem,M.node,q,options{:});
    end
end

numview = getcharin('view',varargin);
up_vector = getcharin('camup',varargin);
camera_position = getcharin('campos',varargin);
if ~isempty(numview)
    view(numview)
elseif dim==3
    view(3)
end
if ~isempty(up_vector)
    camup(up_vector)
end
if ~isempty(camera_position)
    camup(camera_position)
end

if ~(dim==1 && ischarin('courbe',varargin))
    axis off
end

if nargout>=1
    varargout{1} = Helem;
end
