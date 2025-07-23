function varargout = gmshwithfaces(D,varargin)
% function G = gmshwithfaces(D,varargin)
% function u = gmshwithfaces(D,P,varargin)
% D : GEOMOBJECT
% P : POINT

filename = getcharin('filename',varargin,'gmsh_file');
indim = getcharin('indim',varargin,getindim(D));
recombine = ischarin('recombine',varargin);
opencascade = any(ischarin({'occ','OpenCASCADE'},varargin));
varargin = delcharin({'filename','indim'},varargin);
varargin = delonlycharin({'recombine','occ','OpenCASCADE'},varargin);

if D.dim<=1
    if isscalar(varargin{1}) && ~isa(varargin{1},'POINT')
        G = gmshfile(D,varargin{:});
    else
        G = gmshfilewithpoints(D,varargin{:});
    end
else
    if isscalar(varargin{1}) && ~isa(varargin{1},'POINT')
        if opencascade
            [G,number] = gmshfileocc(D,varargin{:});
        else
            [G,number] = gmshfile(D,varargin{:});
        end
    else
        if opencascade
            [G,number] = gmshfilewithpointsocc(D,varargin{:});
        else
            [G,number] = gmshfilewithpoints(D,varargin{:});
        end
    end
    if recombine
        if D.dim==2
            G = recombinesurface(G,number);
        elseif D.dim==3
            G = recombinesurface(G);
        end
    end
    if D.dim==2
        G = createphysicalsurface(G,number,1);
    elseif D.dim==3
        G = createphysicalvolume(G,number,1);
    end
end

G = setfile(G,filename);
n = max(nargout,1);
varargout = cell(1,n);
[varargout{:}] = gmsh2femobject(indim,G,getdim(D):-1:getdim(D)-1);
