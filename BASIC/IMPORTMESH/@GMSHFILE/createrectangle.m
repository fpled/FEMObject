function u = createrectangle(u,point,extents,numbersurface,varargin)
% function u = createrectangle(u,point,extents,numbersurface)
% function u = createrectangle(u,point,extents,numbersurface,radius)

if length(point)~=3
    error('A rectangle is defined by the 3 coordinates of its lower-left corner: point = [x,y,z].')
end
if length(extents)~=2
    error('A rectangle is defined by the 2 extents (width and height): extents = [dx,dy].')
end

u = createentity(u,'Rectangle',[point(:)' extents(:)' varargin{:}],numbersurface);
