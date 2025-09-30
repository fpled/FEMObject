function [u,tag] = extrude(u,vect,name,number,varargin)
% function u = extrude(u,vect,name,number[,tag,'recombine','Layers',layers])
% function [u,tag] = extrude(u,vect,name,number[,tag,'recombine','Layers',layers])
% tag = 'out' by default if nargout==2

vect = vect(:).'; % 3 components of the translation vector
if length(vect)==2, vect = [vect, 0]; end

% Flags for options
recombine = ischarin('recombine',varargin);
layers    = getcharin('Layers',varargin,[]);

varargin = delonlycharin('recombine',varargin);
varargin = delcharin('Layers',varargin);

% Output capture
tag = ''; prefix = '';
if nargin>4 && ~isempty(varargin) && ischar(varargin{1})
    tag = varargin{1}; prefix = [tag '[] = '];
elseif nargout==2
    tag = 'out';       prefix = [tag '[] = '];
end

u = stepcounter(u);

entities = [name tagsintobraces(number) ';'];

options = '';
if ~isempty(layers)
    options = [options ' Layers{' num2str(layers) '};'];
end
if recombine
    options = [options ' Recombine;'];
end

s = [prefix 'Extrude ' tagsintobraces(vect) ...
    ' { ' entities options ' }'];
if ~isempty(prefix), s = [s ';']; end
s = [s '\n'];

u = addstring(u,s);
