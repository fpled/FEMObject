function varargout = adaptmesh(M,q,filename,varargin)
% function varargout = adaptmesh(M,q,filename,varargin)

if ~isa(M,'MODEL')
    varargout = cell(1,nargout);
    [varargout{:}] = adaptmesh(q,M,filename,varargin{:});
    return
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
if verLessThan('matlab','9.1') % compatibility (<R2016b)
    outputFile = ~isempty(strfind(options,'-out'));
else
    outputFile = contains(options,'-out');
end

%% sol file as size map for mesh file
% G = exportfile(G,'.msh','mesh',varargin{:});
% G = writefilesol(G,3,q); % indim=3 for 2D .mesh file generated using Gmsh (coordinates in 3D)
%                          % indim=2 for 2D .mesh file generated using Medit (coordinates in 2D)
% options = [' -sol ' getfilesol(G) ' ' options];
% if ~outputFile
%     options = [options ' -out ' getfile(G,'.mesh')];
% end
% switch dim
%     case 2
%         if verLessThan('matlab','9.1') % compatibility (<R2016b)
%             meshOutput = ~isempty(strfind(options,'.mesh'));
%             opt3dMedit = ~isempty(strfind(options,'-3dMedit'));
%         else
%             meshOutput = contains(options,'.mesh');
%             opt3dMedit = contains(options,'-3dMedit');
%         end
%         if meshOutput && ~opt3dMedit
%             options = [options ' -3dMedit 2']; % to load a 2D .mesh file created with Gmsh (coordinates in 3D) and to produce a Gmsh 2D .mesh file (coordinates in 3D)
%         end
%         G = runfilemmg2d(G,'.mesh',options);
%     case 3
%         G = runfilemmg3d(G,'.mesh',options);
% end
% G = exportfile(G,'.mesh','msh2',varargin{:});

%% NodeData field in msh file
if ~outputFile
    options = [options ' -out ' getfile(G,'.mesh')];
end
if verLessThan('matlab','9.1') % compatibility (<R2016b)
    mshOutput = ~isempty(strfind(options,'.msh'));
else
    mshOutput = contains(options,'.msh');
end
if mshOutput
    G = deletenodedata(G);
end
G = appendnodedata(G,q);
switch dim
    case 2
        if verLessThan('matlab','9.1') % compatibility (<R2016b)
            meshOutput = ~isempty(strfind(options,'.mesh'));
            opt3dMedit = ~isempty(strfind(options,'-3dMedit'));
        else
            meshOutput = contains(options,'.mesh');
            opt3dMedit = contains(options,'-3dMedit');
        end
        if meshOutput && ~opt3dMedit
            options = [options ' -3dMedit 1']; % to force mmg2d to produce a 2D .mesh file readable by Gmsh (not anymore compatible with Medit)
            % options = [options ' -3dMedit 2']; % it also works
        end
        G = runfilemmg2d(G,'.msh',options);
    case 3
        G = runfilemmg3d(G,'.msh',options);
end
G = exportfile(G,'.mesh','msh2',varargin{:});

n=max(nargout,1);
varargout = cell(1,n);
[varargout{:}] = gmsh2femobject(indim,G,dim:-1:dim-n+1,varargin{:});
