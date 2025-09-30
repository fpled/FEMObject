function [u,tag] = affine(u,entries,name,number,varargin)
% function u = affine(u,entries,name,number[,tag,'duplicate'])
% function [u,tag] = affine(u,entries,name,number[,tag,'duplicate'])
% tag = 'out' by default if nargout==2

entries = entries(:).'; % 16 entries of the 4 x 4 affine transformation matrix given by row (only 12 entries can be provided for convenience)

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

s = [prefix 'Affine ' tagsintobraces(entries) ...
    ' { ' entities ' }'];
if ~isempty(prefix), s = [s ';']; end
s = [s '\n'];

u = addstring(u,s);
