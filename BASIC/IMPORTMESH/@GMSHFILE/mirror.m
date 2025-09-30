function [u,tag] = mirror(u,coefficients,name,number,varargin)
% function u = mirror(u,coefficients,name,number[,tag,'duplicate'])
% function [u,tag] = mirror(u,coefficients,name,number[,tag,'duplicate'])
% tag = 'out' by default if nargout==2

coefficients = coefficients(:).'; % 4 coefficients of the plane's equation

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

s = [prefix 'Symmetry ' tagsintobraces(coefficients) ...
    ' { ' entities ' }'];
if ~isempty(prefix), s = [s ';']; end
s = [s '\n'];

u = addstring(u,s);
