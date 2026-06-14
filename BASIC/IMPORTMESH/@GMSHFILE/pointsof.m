function u = pointsof(u,name,number,varargin)
% function u = pointsof(u,name,number[,tag,'duplicate'])

% Duplicate
duplicate = ischarin('duplicate',varargin);
varargin = delonlycharin('duplicate',varargin);

if nargin<3 || isempty(number)
    number = ':';
end

% Output capture
prefix = '';
if nargin>3 && ~isempty(varargin)
    p = tag2str(varargin{1});
    prefix = [p ' = '];
end

u = stepcounter(u);

entities = [name tagsintobraces(number) ';'];
if duplicate
    entities = ['Duplicata{ ' entities ' }'];
end

s = [prefix 'PointsOf{ ' entities ' }'];
if ~isempty(prefix), s = [s ';']; end
s = [s '\n'];

u = addstring(u,s);
