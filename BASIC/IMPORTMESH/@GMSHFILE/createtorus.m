function u = createtorus(u,center,radii,numbervolume,angle)
% function u = createtorus(u,center,radii,numbervolume)
% function u = createtorus(u,center,radii,numbervolume,angle)

if ~(isnumeric(center) && numel(center)==3)
    error('createtorus:CenterInvalid', ...
        'center must be a numeric 1x3 vector [cx,cy,cz].');
end
if ~(isnumeric(radii) && numel(radii)==2)
    error('createtorus:RadiiInvalid', ...
        'radii must be a numeric 1x2 vector [r1,r2].');
end

base = [center(:).' radii(:).'];

if nargin<5 || isempty(angle)
    vals = base;
elseif isnumeric(angle) && isscalar(angle)
    vals = [base angle];
elseif isstring(angle) || ischar(angle)
    vals = [num2cell(base) {tag2str(angle)}];
else
    error('createtorus:AngleInvalid', ...
        'angle must be [], a numeric scalar, or a char/string tag.');
end

u = createentity(u,'Torus',vals,numbervolume);

end
