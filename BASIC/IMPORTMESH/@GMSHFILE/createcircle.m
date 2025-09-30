function u = createcircle(u,center,radius,numbercurve,angles)
% function u = createcircle(u,center,radius,numbercurve)
% function u = createcircle(u,center,radius,numbercurve,angles)

if ~(isnumeric(center) && numel(center)==3)
    error('createcircle:CenterInvalid', ...
        'center must be a numeric 1x3 vector [cx,cy,cz].');
end
if ~(isnumeric(radius) && isscalar(radius))
    error('createcircle:RadiusInvalid', ...
        'radius must be a numeric scalar.');
end

base = [center(:).' radius];

if nargin<5 || isempty(angles)
    vals = base;
    
elseif isnumeric(angles)
    if numel(angles)~=2
        error('createcircle:AnglesNumericSize', ...
            'angles must be a numeric 1x2 vector [angle1,angle2].');
    end
    vals = [base angles(:).'];
    
elseif isstring(angles) && numel(angles)==2
    vals = [num2cell(base) cellstr(angles(:)).'];
    
elseif ischar(angles) && ismatrix(angles) && size(angles,1)==2
    vals = [num2cell(base) cellstr(angles).'];
    
elseif iscell(angles) && numel(angles)==2
    vals = [num2cell(base) cellfun(@tag2str, angles(:).', 'UniformOutput', false)];
    
else
    error('createcircle:AnglesFormatUnsupported', ...
        'Unsupported "angles" format. Use numeric [a1,a2], string array ["a1","a2"], 2xN char matrix [''a1'';''a2''], cell array of character vectors {''a1'',''a2''}, or cell array of numeric/char/string tags {a1,a2}.');
end

u = createentity(u,'Circle',vals,numbercurve);

end
