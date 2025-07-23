function u = createcylinder(u,center,axis,radius,numbervolume,varargin)
% function u = createcylinder(u,center,axis,radius,numbervolume)
% function u = createcylinder(u,center,axis,radius,numbervolume,angle)

if length(center)~=3
    error('A cylinder is defined by the 3 coordinates of its base center: center = [cx,cy,cz].')
end
if isempty(radius) || ~isscalar(radius)
    error('A cylinder is defined by its scalar radius.');
end
if length(axis)~=3
    error('A cylinder is defined by the 3 components of its axis vector: axis = [dx,dy,dz].')
end

u = createentity(u,'Cylinder',[center(:)' axis(:)' radius varargin{:}],numbervolume);
