function u = createcone(u,center,axis,radii,numbervolume,varargin)
% function u = createcone(u,center,axis,radii,numbervolume)
% function u = createcone(u,center,axis,radii,numbervolume,angle)

if length(center)~=3
    error('A cone is defined by the 3 coordinates of its base center: center = [cx,cy,cz]')
end
if length(axis)~=3
    error('A cone is defined by the 3 components of its axis vector: axis = [dx,dy,dz]')
end
if length(radii)~=2
    error('A cone is defined by the 2 radii of its faces: radii = [r1,r2]')
end

u = createentity(u,'Cone',[center(:)' axis(:)' radii(:)' varargin{:}],numbervolume);
