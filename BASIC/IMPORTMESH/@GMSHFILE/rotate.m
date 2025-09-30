function [u,tag] = rotate(u,direction,point,angle,name,number,varargin)
% function u = rotate(u,direction,point,angle,name,number[,tag,'duplicate'])
% function [u,tag] = rotate(u,axis,point,angle,name,number[,tag,'duplicate'])
% tag = 'out' by default if nargout==2

direction = direction(:).'; % 3 components of the direction vector of the rotation axis
point = point(:).';         % 3 coordinates of any point on the rotation axis
if length(direction)==2,  direction  = [direction, 0];  end
if length(point)==2,      point      = [point, 0];      end

% Duplicate
duplicate = ischarin('duplicate',varargin);
varargin = delonlycharin('duplicate',varargin);

% Output capture
tag = ''; prefix = '';
if nargin>6 && ~isempty(varargin) && ischar(varargin{1})
    tag = varargin{1}; prefix = [tag '[] = '];
elseif nargout==2
    tag = 'out';       prefix = [tag '[] = '];
end

u = stepcounter(u);

entities = [name tagsintobraces(number) ';'];
if duplicate
    entities = ['Duplicata{ ' entities ' }'];
end

s = [prefix 'Rotate { ' tagsintobraces(direction) ', ' tagsintobraces(point) ', ' tag2str(angle) ' }' ...
    ' { ' entities ' }'];
if ~isempty(prefix), s = [s ';']; end
s = [s '\n'];

u = addstring(u,s);
