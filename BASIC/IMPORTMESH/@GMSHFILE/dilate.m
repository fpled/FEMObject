function [u,tag] = dilate(u,center,factor,name,number,varargin)
% function u = dilate(u,center,factor,name,number[,tag,'duplicate'])
% function [u,tag] = dilate(u,center,factor,name,number[,tag,'duplicate'])
% tag = 'out' by default if nargout==2

center = center(:).'; % 3 coordinates of the homothety center (center of the homothetic transformation)
if isscalar(factor)
    factor = num2str(factor); % 1 isotropic scaling factor
else
    factor = factor(:).';     % 3 anisotropic scaling factors along the 3 coordinate axes
    factor = tagsintobraces(factor);
end

% Duplicate
duplicate = ischarin('duplicate',varargin);
varargin = delonlycharin('duplicate',varargin);

% Output capture
tag = ''; prefix = '';
if nargin>5 && ~isempty(varargin) && ischar(varargin{1})
    tag = varargin{1}; prefix = [tag '[] = '];
elseif nargout==2
    tag = 'out';       prefix = [tag '[] = '];
end

u = stepcounter(u);

entities = [name tagsintobraces(number) ';'];
if duplicate
    entities = ['Duplicata{ ' entities ' }'];
end

s = [prefix 'Dilate { ' tagsintobraces(center) ', ' factor ' }' ...
    ' { ' entities ' }'];
if ~isempty(prefix), s = [s ';']; end
s = [s '\n'];

u = addstring(u,s);
