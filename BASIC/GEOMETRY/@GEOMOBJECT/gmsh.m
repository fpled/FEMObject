function varargout = gmsh(D,varargin)
% function u = gmsh(D,varargin)
% function u = gmsh(D,P,varargin)
% D : GEOMOBJECT
% P : POINT

filename = getcharin('filename',varargin,'gmsh_file');
indim = getcharin('indim',varargin,getindim(D));
recombine = ischarin('recombine',varargin);
varargin = delcharin({'filename','indim'},varargin);
varargin = delonlycharin('recombine',varargin);

if D.dim<=1
    if isscalar(varargin{1}) && ~isa(varargin{1},'POINT')
        G = gmshfile(D,varargin{:});
    else
        G = gmshfilewithpoints(D,varargin{:});
    end
else
    if isscalar(varargin{1}) && ~isa(varargin{1},'POINT')
        [G,numbersurface] = gmshfile(D,varargin{:});
    else
        [G,numbersurface] = gmshfilewithpoints(D,varargin{:});
    end
    if recombine
        if D.dim==2
            G = recombinesurface(G,numbersurface);
        elseif D.dim==3
            G = recombinesurface(G);
        end
    end
    G = createphysicalsurface(G,numbersurface,1);
    if D.dim==3
        G = createphysicalvolume(G,1,1);
    end
end

G = setfile(G,filename);
n=max(nargout,1);
varargout=cell(1,n);
[varargout{:}] = gmsh2femobject(indim,G,getdim(D):-1:getdim(D)-n+1);    
 