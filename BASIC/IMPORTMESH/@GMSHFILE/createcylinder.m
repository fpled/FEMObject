function u = createcylinder(u,center,axis,radius,numbervolume,angle)
% function u = createcylinder(u,center,axis,radius,numbervolume)
% function u = createcylinder(u,center,axis,radius,numbervolume,angle)

if ~(isnumeric(center) && numel(center)==3)
    error('createcylinder:CenterInvalid', ...
        'center must be a numeric 1x3 vector [cx,cy,cz].');
end
if ~(isnumeric(axis) && numel(axis)==3)
    error('createcylinder:AxisInvalid', ...
        'axis must be a numeric 1x3 vector [dx,dy,dz].');
end
if ~(isnumeric(radius) && isscalar(radius))
    error('createcylinder:RadiusInvalid', ...
        'radius must be a numeric scalar.');
end

base = [center(:).' axis(:).' radius];

if nargin<6 || isempty(angle)
    vals = base;
elseif isnumeric(angle) && isscalar(angle)
    vals = [base angle];
elseif isstring(angle) || ischar(angle)
    vals = [num2cell(base) {tag2str(angle)}];
else
    error('createcylinder:AngleInvalid', ...
        'angle must be [], a numeric scalar, or a char/string tag.');
end

u = createentity(u,'Cylinder',vals,numbervolume);

end
