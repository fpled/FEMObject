function varargout = adaptmesh(M,q,filename,varargin)
% function varargout = adaptmesh(M,q,filename,varargin)

if ~isa(M,'MODEL')
    varargout = cell(1,nargout);
    [varargout{:}] = adaptmesh(q,M,filename,varargin{:});
    return
end

if verLessThan('matlab','9.1') % compatibility (<R2016b)
    contain = @(str,pat) ~isempty(strfind(str,pat));
else
    contain = @contains;
end

if israndom(q)
    error('la solution est aleatoire')
end
q = unfreevector(M,q);

indim = getindim(M);
dim = getdim(M);

G = GMSHFILE();
if nargin>=3 && ischar(filename)
    G = setfile(G,filename);
end

options = getcharin('mmgoptions',varargin);

%% sol file as size map for mesh file
% G = exportfile(G,'.msh','mesh',varargin{:});
% G = writefilesol(G,3,q); % indim=3 for 2D .mesh file generated using Gmsh (coordinates in 3D), indim=2 for 2D .mesh file generated using Medit (coordinates in 2D)
% options = [' -sol ' getfilesol(G) ' ' options];
% if ~contain(options,'-out')
%     % options = [options ' -out ' getfile(G,'.mesh')]; % output .mesh file
%     options = [options ' -out ' getfile(G,'.msh')]; % output .msh file
% end
% 
% switch dim
%     case 2
%         if ~contain(options,'3dMedit')
%             options = [options ' -3dMedit 2']; % to load a 2D .mesh file created with Gmsh (coordinates in 3D) and to produce a Gmsh 2D .mesh file (coordinates in 3D)
%         end
%         G = runfilemmg2d(G,'.mesh',options);
%     case 3
%         G = runfilemmg3d(G,'.mesh',options);
% end
% 
% if contain(options,'.mesh')
%     G = exportfile(G,'.mesh','msh2',varargin{:});
% elseif contain(options,'.msh')
%     G = exportfile(G,'.msh','msh2',varargin{:});
%     G = exportfile(G,'.msh','msh2',varargin{:}); % export twice for renumbering of mesh nodes/elements
% end

%% NodeData field in msh file
if ~contain(options,'-out')
    % options = [options ' -out ' getfile(G,'.mesh')]; % output .mesh file
    options = [options ' -out ' getfile(G,'.msh')]; % output .msh file
end
G = updatenodedata(G,q);
switch dim
    case 2
        if contain(options,'.mesh') && ~contain(options,'-3dMedit')
            options = [options ' -3dMedit 1']; % to force mmg2d to produce a 2D .mesh file readable by Gmsh (not anymore compatible with Medit)
            % options = [options ' -3dMedit 2']; % it also works
        end
        G = runfilemmg2d(G,'.msh',options);
    case 3
        G = runfilemmg3d(G,'.msh',options);
end
if contain(options,'.mesh')
    G = exportfile(G,'.mesh','msh2',varargin{:});
end

n=max(nargout,1);
varargout = cell(1,n);
[varargout{:}] = gmsh2femobject(indim,G,dim:-1:dim-n+1,varargin{:});
