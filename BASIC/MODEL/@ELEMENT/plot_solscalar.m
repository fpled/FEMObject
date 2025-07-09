function varargout = plot_solscalar(elem,node,q,varargin)
% function varargout = plot_solscalar(elem,node,q,varargin)

nodecoord = double(getcoord(node));
dim = size(nodecoord,2);
connec = calc_conneclocal(elem,node) ;
qe = localize(elem,q);
xnode = node(elem);
ampl = getcharin('ampl',varargin,0);
varargin = delcharin('ampl',varargin);
numview = getcharin('view',varargin);
up_vector = getcharin('camup',varargin);
camera_position = getcharin('campos',varargin);

if ~isenrich(elem)
    u = reshape(double(qe),[getnbnode(elem),getnbelem(elem)]);
else
    p = permute(nodelocalcoord(elem),[4,2,3,1]);
    u = double(calc_N(elem,xnode,p)*qe);
    u = permute(u,[1,2,4,3]);
end

globconnec = zeros(size(connec))';
globconnec(:) = 1:numel(connec);
globconnec = globconnec';
% globconnec = reshape(1:numel(connec),size(connec));
faces = patchfaces(elem,globconnec);
nbfaces = size(patchfaces(elem,globconnec(1,:)),1);

connec = connec';
vertpos = nodecoord(connec(:),:);

if ischarin('sigma',varargin) || ischarin('epsilon',varargin) || ischarin('energyint',varargin)
    mat = getmaterial(elem);
    mid = nodelocalcoord(elem);
    mid = sum(mid,1)/size(mid,1);
    if ischarin('sigma',varargin)
        se = sigma(mat,elem,xnode,mid,qe);
        k = getcharin('sigma',varargin);
    elseif ischarin('epsilon',varargin)
        se = calc_B(elem,xnode,mid)*qe;
        k = getcharin('epsilon',varargin);
    elseif ischarin('energyint',varargin)
        k = getcharin('energyint',varargin);
        se = energyint(mat,elem,xnode,mid,qe,k);
        k = 1;
    end
    se = double(sigmacompo(se,k,elem));
    se = se(:);
    se = repmat(se',nbfaces,1);
    varargin = delcharin('sigma',varargin);
    varargin = delcharin('epsilon',varargin);
    varargin = delcharin('energyint',varargin);
    varargin = setcharin('facecolor',varargin,'flat');
    varargin = setcharin('facevertexcdata',varargin,se(:));
    
    if ischarin('surface',varargin)
        zpos = repmat(se(:)',getnbnode(elem),1);
        vertpos = [vertpos,zpos(:)];
        varargin = delonlycharin('surface',varargin);
    end
    
    H = patch('faces',faces,'vertices',vertpos,varargin{:});
    axis image
else
    if strcmp(getcharin('facecolor',varargin),'none')
        varargin = setcharin('facecolor',varargin,'interp');
    end
    if strcmp(getcharin('edgecolor',varargin),'none')
        varargin = setcharin('edgecolor',varargin,'interp');
    end
    varargin = delcharin('displ',varargin);
    varargin = delcharin('rotation',varargin);
    
    if (ischarin('surface',varargin) || ischarin('surfacemesh',varargin))% && getdim(elem)==2
        varargin = delonlycharin('surface',varargin);
        varargin = setcharin('facevertexcdata',varargin,q(:));
        H = surface(elem,node,'facevertexcdata',q,varargin{:});
        
        if ischarin('surfacemesh',varargin)
            varargin = setcharin('edgecolor',varargin,'k');
            varargin = setcharin('facecolor',varargin,'none');
            varargin = delonlycharin('surfacemesh',varargin);
            plot(elem,node,'facevertexcdata',q,varargin{:});
        end
        axis image
        axis square
        if ~isempty(numview)
            view(numview)
        else
            view(3)
        end
    else
        varargin = setcharin('facevertexcdata',varargin,u(:));
        H = patch('faces',faces,'vertices',vertpos,varargin{:});
        axis image
    end
end

if ~isempty(up_vector)
    camup(up_vector)
end
if ~isempty(camera_position)
    campos(camera_position);
end

if nargout>=1
    varargout{1} = H;
end


function se = sigmacompo(se,ksigma,elem)

if isa(ksigma,'char') && strcmp(ksigma,'mises')
    switch getindim(elem)
        case 1
            se = se(1);
        case 2
            tracese = 1/3*(se(1)+se(2));
            se(1) = se(1) - tracese;
            se(2) = se(2) - tracese;
            se = sqrt(3/2*(se(1)*se(1) + se(2)*se(2) + 2*se(3)*se(3)));
        case 3
            tracese = 1/3*(se(1)+se(2)+se(3));
            se(1) = se(1) - tracese;
            se(2) = se(2) - tracese;
            se(3) = se(3) - tracese;
            se = sqrt(3/2*(se(1)*se(1) + se(2)*se(2) + se(3)*se(3)...
                + 2*(se(4)*se(4) + se(5)*se(5) + se(6)*se(6))));
    end
else
    se = se(ksigma);
end

return
