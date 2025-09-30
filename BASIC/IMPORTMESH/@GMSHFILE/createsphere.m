function u = createsphere(u,center,radius,numbervolume,varargin)
% function u = createsphere(u,center,radius,numbervolume)
% function u = createsphere(u,center,radius,numbervolume,angles)
% function u = createsphere(u,center,radius,numbervolume,angle1,angle2,angle3)
% 
% For partial spheres:
%   angles = [angle1, angle2, angle3]
%   angle1, angle2 : polar opening in [-pi/2, +pi/2]
%   angle3         : azimuthal opening in [0, 2*pi]

if ~(isnumeric(center) && numel(center)==3)
    error('createsphere:CenterInvalid', ...
        'center must be a numeric 1x3 vector [cx,cy,cz].');
end
if ~(isnumeric(radius) && isscalar(radius))
    error('createsphere:RadiusInvalid', ...
        'radius must be a numeric scalar.');
end
if numel(varargin)>3
    error('createsphere:TooManyAngles', ...
        'At most 1 vector "angles" or 3 scalars "angle1,angle2,angle3" are allowed.');
end

base = [center(:).' radius];

if isempty(varargin)
    vals = base;
    
elseif isscalar(varargin)
    angles = varargin{1};
    
    if isnumeric(angles)
        if numel(angles)~=3
            error('createsphere:AnglesNumericSize', ...
                'angles must be a numeric 1x3 vector [angle1,angle2,angle3].');
        end
        vals = [base angles(:).'];
        
    elseif isstring(angles) && numel(angles)==3
        vals = [num2cell(base) cellstr(angles(:)).'];
        
    elseif ischar(angles) && ismatrix(angles) && size(angles,1)==3
        vals = [num2cell(base) cellstr(angles).'];
        
    elseif iscell(angles) && numel(angles)==3
        vals = [num2cell(base) cellfun(@tag2str, angles(:).', 'UniformOutput', false)];
        
    else
        error('createsphere:AnglesFormatUnsupported', ...
            'Unsupported "angles" format. Use numeric [a1,a2,a3], string array ["a1","a2","a3"], 3xN char matrix [''a1'';''a2'';''a3''], cell array of character vectors {''a1'',''a2'',''a3''}, or cell array of numeric/char/string tags {a1,a2,a3}.');
    end
    
else
    if all(cellfun(@(x) isnumeric(x) && isscalar(x), varargin))
        vals = [base cell2mat(varargin)];
    else
        vals = [num2cell(base) cellfun(@tag2str, varargin, 'UniformOutput', false)];
    end
end

u = createentity(u,'Sphere',vals,numbervolume);

end
