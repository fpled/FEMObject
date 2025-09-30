function [u,tag] = translate(u,vect,name,number,varargin)
% function u = translate(u,vect,name,number[,tag,'duplicate'])
% function [u,tag] = translate(u,vect,name,number[,tag,'duplicate'])
% tag = 'out' by default if nargout==2

vect = vect(:).'; % 3 components of the translation vector
if length(vect)==2, vect = [vect, 0]; end

% Duplicate
duplicate = ischarin('duplicate',varargin);
varargin = delonlycharin('duplicate',varargin);

% Output capture
tag = ''; prefix = '';
if nargin>4 && ~isempty(varargin) && ischar(varargin{1})
    tag = varargin{1}; prefix = [tag '[] = '];
elseif nargout==2
    tag = 'out';       prefix = [tag '[] = '];
end

u = stepcounter(u);

entities = [name tagsintobraces(number) ';'];
if duplicate
    entities = ['Duplicata{ ' entities ' }'];
end

s = [prefix 'Translate ' tagsintobraces(vect) ...
    ' { ' entities ' }'];
if ~isempty(prefix), s = [s ';']; end
s = [s '\n'];

u = addstring(u,s);
