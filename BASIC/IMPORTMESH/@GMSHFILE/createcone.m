function u = createcone(u,center,axis,radii,numbervolume,angle)
% function u = createcone(u,center,axis,radii,numbervolume)
% function u = createcone(u,center,axis,radii,numbervolume,angle)

if ~(isnumeric(center) && numel(center)==3)
    error('createcone:CenterInvalid', ...
        'center must be a numeric 1x3 vector [cx,cy,cz]');
end
if ~(isnumeric(axis) && numel(axis)==3)
    error('createcone:AxisInvalid', ...
        'axis must be a numeric 1x3 vector [dx,dy,dz]');
end
if ~(isnumeric(radii) && numel(radii)==2)
    error('createcone:RadiiInvalid', ...
        'radii must be a numeric 1x2 vector [r1,r2]');
end

base = [center(:).' axis(:).' radii(:).'];

if nargin<6 || isempty(angle)
    vals = base;
elseif isnumeric(angle) && isscalar(angle)
    vals = [base angle];
elseif isstring(angle) || ischar(angle)
    vals = [num2cell(base) {tag2str(angle)}];
else
    error('createcone:AngleInvalid', ...
        'angle must be [], a numeric scalar, or a char/string tag.');
end

u = createentity(u,'Cone',vals,numbervolume);

end
