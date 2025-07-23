function u = createsphere(u,center,radius,numbervolume,varargin)
% function u = createsphere(u,center,radius,numbervolume)
% function u = createsphere(u,center,radius,numbervolume,angle1,angle2,angle3)

if length(center)~=3
    error('A sphere is defined by the 3 coordinates of its center: center = [x,y,z].')
end
if isempty(radius) || ~isscalar(radius)
    error('A sphere is defined by its scalar radius.');
end

u = createentity(u,'Sphere',[center(:)' radius varargin{:}],numbervolume);
