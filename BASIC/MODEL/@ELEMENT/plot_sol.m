function varargout = plot_sol(elem,node,q,varargin)
% function varargout = plot_sol(elem,node,q,varargin)

mat = getmaterial(elem);
if ~isa(mat,'ELAS_ISOT') && ~isa(mat,'ELAS_ISOT_TRANS') && ~isa(mat,'ELAS_ANISOT') && ~isa(mat,'ELAS_BEAM') && ~isa(mat,'ELAS_BEAM_ISOT_TRANS') && ~isa(mat,'ELAS_SHELL') && ~isa(mat,'ELAS_SHELL_ISOT_TRANS') && getnbddlpernode(elem)==1
    varargout = cell(1,nargout);
    [varargout{:}] = plot_solscalar(elem,node,q,varargin{:});
    return
end

nodecoord = double(getcoord(node));
dim = size(nodecoord,2);
connec = calc_conneclocal(elem,node);
qe = localize(elem,q);
xnode = node(elem);
ampl = getcharin('ampl',varargin,0);
varargin = delcharin('ampl',varargin);

if ~isenrich(elem)
    u = reshape(double(qe),[getnbddlpernode(elem),getnbnode(elem),getnbelem(elem)]);
    u = permute(u,[1,4,2,3]);
else
    p = permute(nodelocalcoord(elem),[4,2,3,1]);
    u = double(calc_N(elem,xnode,p)*qe);
    u = permute(u,[1,2,4,3]);
end
u = reshape(u,[size(u,1),numel(u)/size(u,1)])';

globconnec = zeros(size(connec))';
globconnec(:) = 1:numel(connec);
globconnec = globconnec';
% globconnec = reshape(1:numel(connec),size(connec));
faces = patchfaces(elem,globconnec);
nbfaces = size(patchfaces(elem,globconnec(1,:)),1);

connec = connec';
if ~isa(mat,'ELAS_BEAM') && ~isa(mat,'ELAS_BEAM_ISOT_TRANS') && ~isa(mat,'ELAS_SHELL') && ~isa(mat,'ELAS_SHELL_ISOT_TRANS')
    vertpos = nodecoord(connec(:),:)+ampl*u;
else
    vertpos = nodecoord(connec(:),:)+ampl*u(:,1:dim);
end

if ischarin('sigma',varargin) || ischarin('epsilon',varargin) || ischarin('energyint',varargin)
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
    varargin = setcharin('FaceColor',varargin,'flat');
    varargin = setcharin('FaceVertexCData',varargin,se(:));
    
    if ischarin('surface',varargin)
        zpos = repmat(se(:)',getnbnode(elem),1);
        vertpos = [vertpos,zpos(:)];
        varargin = delonlycharin('surface',varargin);
    end
else
    if strcmp(getcharin('FaceColor',varargin),'none')
        varargin = setcharin('FaceColor',varargin,'interp');
    end
    if (isa(mat,'ELAS_BEAM') || isa(mat,'ELAS_BEAM_ISOT_TRANS')) && strcmp(getcharin('EdgeColor',varargin),'none')
        varargin = setcharin('EdgeColor',varargin,'interp');
    end
    if ischarin('displ',varargin)
        k = getcharin('displ',varargin,1);
    elseif ischarin('rotation',varargin)
        k = dim + getcharin('rotation',varargin,1);
    else
        k = 1;
    end
    varargin = delcharin('displ',varargin);
    varargin = delcharin('rotation',varargin);
    varargin = setcharin('FaceVertexCData',varargin,u(:,k));
end

H = patch('Faces',faces,'Vertices',vertpos,varargin{:});

axis image

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
