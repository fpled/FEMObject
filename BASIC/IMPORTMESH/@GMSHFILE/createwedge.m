function u = createwedge(u,point,extents,numbervolume,varargin)
% function u = createwedge(u,point,extents,numbervolume)
% function u = createwedge(u,point,extents,numbervolume,topXextent)

if length(center)~=3
    error('A wedge is defined by the 3 coordinates of its right-angle point: point = [x,y,z].')
end
if length(extents)~=3
    error('A wedge is defined by the 3 extents (width, height and depth): extents = [dx,dy,dz].')
end

u = createentity(u,'Wedge',[point(:)' extents(:)' varargin{:}],numbervolume);
