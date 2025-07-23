function u = createbox(u,point,extents,numbervolume)
% function u = createbox(u,point,extents,numbervolume)

if length(point)~=3
    error('A box is defined by the 3 coordinates of its lower-left corner: point = [x,y,z].')
end
if length(extents)~=3
    error('A box is defined by the 3 extents (width, height and depth): extents = [dx,dy,dz].')
end

u = createentity(u,'Box',[point(:)' extents(:)'],numbervolume);
