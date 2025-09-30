function u = createwedge(u,point,extents,numbervolume,ltx)
% function u = createwedge(u,point,extents,numbervolume)
% function u = createwedge(u,point,extents,numbervolume,ltx)

if ~(isnumeric(point) && numel(point)==3)
    error('createwedge:PointInvalid', ...
        'point must be a numeric 1x3 vector [x,y,z].');
end
if ~(isnumeric(extents) && numel(extents)==3)
    error('createwedge:ExtentsInvalid', ...
        'extents must be a numeric 1x3 vector [dx,dy,dz].');
end

base = [point(:).' extents(:).'];

if nargin<5 || isempty(ltx)
    vals = base;
else
    if ~(isnumeric(ltx) && isscalar(ltx))
        error('createwedge:LtxInvalid', ...
            'ltx must be a numeric scalar.');
    end
    vals = [base ltx];
end

u = createentity(u,'Wedge',vals,numbervolume);

end
