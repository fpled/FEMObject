function u = createbox(u,point,extents,numbervolume)
% function u = createbox(u,point,extents,numbervolume)

if ~(isnumeric(point) && numel(point)==3)
    error('createbox:PointInvalid', ...
        'point must be a numeric 1x3 vector [x,y,z].');
end
if ~(isnumeric(extents) && numel(extents)==3)
    error('createbox:ExtentsInvalid', ...
        'extents must be a numeric 1x3 vector [dx,dy,dz].');
end

vals = [point(:).' extents(:).'];

u = createentity(u,'Box',vals,numbervolume);

end
