function [u,tag] = revolve(u,direction,point,angle,name,number,varargin)
% function u = revolve(u,direction,point,angle,name,number[,tag,'recombine','Layers',layers])
% function [u,tag] = revolve(u,direction,point,angle,name,number[,tag,'recombine','Layers',layers])
% tag = 'out' by default if nargout==2

direction = direction(:).'; % 3 components of the direction vector of the rotation axis
point = point(:).';         % 3 coordinates of any point on the rotation axis
if length(direction)==2, direction = [direction, 0]; end
if length(point)==2,     point     = [point, 0];     end
if isempty(angle), angle = '2*Pi'; end

% Flags for options
recombine = ischarin('recombine',varargin);
layers    = getcharin('Layers',varargin,[]);

varargin = delonlycharin('recombine',varargin);
varargin = delcharin('Layers',varargin);

% Output capture
tag = ''; prefix = '';
if nargin>6 && ~isempty(varargin) && ischar(varargin{1})
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

s = [prefix 'Extrude { ' tagsintobraces(direction) ', ' tagsintobraces(point) ', ' tag2str(angle) ' }' ...
    ' { ' entities options ' }'];
if ~isempty(prefix), s = [s ';']; end
s = [s '\n'];

u = addstring(u,s);
