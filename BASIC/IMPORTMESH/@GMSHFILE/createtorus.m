function u = createtorus(u,center,radii,numbervolume,varargin)
% function u = createtorus(u,center,radii,numbervolume)
% function u = createtorus(u,center,radii,numbervolume,angle)

if length(center)~=3
    error('A torus is defined by the 3 coordinates of its center: center = [cx,cy,cz].')
end
if length(radii)~=2
    error('A torus is defined by its 2 radii: radii = [r1,r2].')
end

u = createentity(u,'Torus',[center(:)' radii(:)' varargin{:}],numbervolume);
