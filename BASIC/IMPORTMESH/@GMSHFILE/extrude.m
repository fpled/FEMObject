function [u,out] = extrude(u,vect,name,number,varargin)
% function u = extrude(u,vect,name,number)
% function [u,out] = extrude(u,vect,name,number)

u = stepcounter(u);
vect = vect(:)';
if length(vect)==2
    vect = [vect,0];
end

% Handle 'Recombine' option
if ischarin('recombine',varargin)
    recombine = ' Recombine;';
else
    recombine = '';
end
varargin = delonlycharin('recombine',varargin);

% Handle output variable
out = '';
prefix = '';
if nargin>4 && ~isempty(varargin) && ischar(varargin{1})
    out = varargin{1};
    prefix = [out '[] = '];
elseif nargout==2
    out = 'out';
    prefix = [out '[] = '];
end

s = [prefix 'Extrude' valuesintobraces(vect) ...
    ' { ' name '{' num2str(number) '} ;' recombine ' }'];
if ~isempty(prefix)
    s = [s ';'];
end
s = [s '\n'];

u = addstring(u,s);
