function u = createrectangle(u,point,extents,numbersurface,roundedradius)
% function u = createrectangle(u,point,extents,numbersurface)
% function u = createrectangle(u,point,extents,numbersurface,roundedradius)

if ~(isnumeric(point) && numel(point)==3)
    error('createrectangle:PointInvalid', ...
        'point must be a numeric 1x3 vector [x,y,z].');
end
if ~(isnumeric(extents) && numel(extents)==2)
    error('createrectangle:ExtentsInvalid', ...
        'extents must be a numeric 1x2 vector [dx,dy].');
end

base = [point(:).' extents(:).'];

if nargin<5 || isempty(roundedradius)
    vals = base;
else
    if ~(isnumeric(roundedradius) && isscalar(roundedradius))
        error('createrectangle:RoundedRadiusInvalid', ...
            'roundedradius must be a numeric scalar.');
    end
    vals = [base roundedradius];
end

u = createentity(u,'Rectangle',vals,numbersurface);

end