function f = getfilesol(u,suffix)
% function u = getfilesol(u,suffix)

if nargin==1 || isempty(suffix)
    f = getfile(u,'.sol');
else
    f = getfile(u,'.sol',suffix);
end
