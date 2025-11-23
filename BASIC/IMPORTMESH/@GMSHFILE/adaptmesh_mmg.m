function u = adaptmesh_mmg(u,q,dim,varargin)
% function u = adaptmesh_mmg(u,q,dim)
% mesh adaptation with MMG software
% 1) using a .sol file as size map for .mesh file (Medit file format)
% export u.msh to u.mesh
% write u.sol
% launch mmg2d_O3 u.mesh -sol u.sol -out u.mesh -3dMedit 2 % for dim=2
% launch mmg3d_O3 u.mesh -sol u.sol -out u.mesh            % for dim=3
% export u.mesh to u.msh
% 2) using a NodeData field in .msh file (Gmsh file format)
% append or update NodeData field in .msh file
% launch mmg2d_O3 u.msh -out u.msh % for dim=2
% launch mmg3d_O3 u.msh -out u.msh % for dim=3

% function u = adaptmesh_mmg(u,q,dim,'mmgoptions',mmgoptions,'gmshoptions',gmshoptions)
% mmgoptions: char containing mmg options
% gmshoptions: char containing gmsh options

if nargin==3
    error('rentrer les options de mmg')
end

if verLessThan('matlab','9.1') % compatibility (<R2016b)
    contain = @(str,pat) ~isempty(strfind(str,pat));
else
    contain = @contains;
end

options = getcharin('mmgoptions',varargin);

%% Use a .sol file as size map for .mesh file (Medit file format)
% u = exportfile(u,'.msh','mesh',varargin{:});
% 
% u = writefilesol(u,3,q); % indim=3 for 2D .mesh file generated using Gmsh (coordinates in 3D), indim=2 for 2D .mesh file generated using Medit (coordinates in 2D)
% 
% options = [' -sol ' getfile(u,'.sol') ' ' options];
% 
% if ~contain(options,'-out')
%     options = [options ' -out ' getfile(u,'.mesh')]; % output .mesh file
%     % options = [options ' -out ' getfile(u,'.msh')]; % output .msh file
% end
% 
% switch dim
%     case 2
%         if ~contain(options,'3dMedit')
%             options = [options ' -3dMedit 2']; % to load a 2D .mesh file created with Gmsh (coordinates in 3D) and to produce a Gmsh 2D .mesh file (coordinates in 3D)
%         end
%         u = runfilemmg2d(u,'.mesh',options);
%     case 3
%         u = runfilemmg3d(u,'.mesh',options);
% end
% 
% if contain(options,'.mesh')
%     u = exportfile(u,'.mesh','msh2',varargin{:});
% elseif contain(options,'.msh')
%     u = exportfile(u,'.msh','msh2',varargin{:});
%     u = exportfile(u,'.msh','msh2',varargin{:}); % export twice for renumbering of mesh nodes/elements
% end

%% Use a NodeData field in .msh file (gmsh file format)
u = updatenodedata(u,q);

if ~contain(options,'-out')
    % options = [options ' -out ' getfile(u,'.mesh')]; % output .mesh file
    options = [options ' -out ' getfile(u,'.msh')]; % output .msh file
end

switch dim
    case 2
        if contain(options,'.mesh') && ~contain(options,'-3dMedit')
            options = [options ' -3dMedit 1']; % to force mmg2d to produce a 2D .mesh file readable by Gmsh (not anymore compatible with Medit)
            % options = [options ' -3dMedit 2']; % it also works
        end
        u = runfilemmg2d(u,'.msh',options);
    case 3
        u = runfilemmg3d(u,'.msh',options);
end

if contain(options,'.mesh')
    u = exportfile(u,'.mesh','msh2',varargin{:});
elseif contain(options,'.msh') && ischarin('export',varargin)
    u = exportfile(u,'.msh','msh2',varargin{:});
end
