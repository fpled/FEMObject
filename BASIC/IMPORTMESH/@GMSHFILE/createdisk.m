function u = createdisk(u,center,radius,numbersurface,radiusy)
% function u = createdisk(u,center,radius,numbersurface)
% function u = createdisk(u,center,radius,numbersurface,radiusy)

if ~(isnumeric(center) && numel(center)==3)
    error('createdisk:CenterInvalid', ...
        'center must be a numeric 1x3 vector [cx,cy,cz].');
end
if ~(isnumeric(radius) && isscalar(radius))
    error('createdisk:RadiusInvalid', ...
        'radius must be a numeric scalar.');
end

base = [center(:).' radius];

if nargin<5 || isempty(radiusy)
    vals = base;
else
    if ~(isnumeric(radiusy) && isscalar(radiusy))
        error('createdisk:RadiusYInvalid', ...
            'radiusy must be a numeric scalar.');
    end
    vals = [base radiusy];
end

u = createentity(u,'Disk',vals,numbersurface);

end
