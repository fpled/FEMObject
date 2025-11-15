function u = adaptmesh_gmsh(u,q,options,varargin)
% function u = adaptmesh_gmsh(u,q,dim)
% launch gmsh u.geo -dim -bgm u.msh

% function u = adaptmesh(u,q,dim,'gmshoptions',options)
% options: char containing gmsh options
% launch gmsh u.geo -dim -bgm u.msh options

if nargin==2
    error('rentrer les options de gsmh')
end

if verLessThan('matlab','9.1') % compatibility (<R2016b)
    contain = @(str,pat) ~isempty(strfind(str,pat));
else
    contain = @contains;
end

if ~isempty(options) && isa(options,'double')
    options = [' -' num2str(options)];
else
    options = '';
end

if ischarin('gmshoptions',varargin)
    options = [options ' ' getcharin('gmshoptions',varargin)];
end

if ~contain(options,'-format')
    options = [options ' -format msh2'];
end

% Load the post-processing view (NodeData field) contained in the .msh file
% (size-map mesh) and apply it as the current background mesh
if ~contain(options,'-bgm')
    options = [options ' -bgm ' getfile(u,'.msh')];
end

% Optional: Use ONLY the background mesh, and disregard any other size constraints (point sizes, curvature, etc.), to define mesh sizes
% if ~contains(options,'-setnumber')
%     options = [options ' -setnumber Mesh.MeshSizeFromPoints 0'];
%     options = [options ' -setnumber Mesh.MeshSizeFromCurvature 0'];
%     options = [options ' -setnumber Mesh.MeshSizeExtendFromBoundary 0'];
% end

% Optional: For 2D mesh, use the "Delaunay" algorithm (Mesh.Algorithm = 5)
% instead of the default "Frontal-Delaunay" algorithm (Mesh.Algorithm = 6)
% While the default "Frontal-Delaunay" 2D meshing algorithm (Mesh.Algorithm = 6)
% usually leads to the highest quality meshes, the "Delaunay" algorithm (Mesh.Algorithm = 5)
% tends to behave better for complex mesh size fields - in particular size fields with large element size gradients.
% if ~contains(options,'-algo') && dim==2
%     options = [options ' -algo del2d'];
% end

u = updatenodedata(u,q);

u = runfile(u,'.geo',options);
