function f = getfilemesh(u,suffix)
% function u = getfilemesh(u,suffix)

if nargin==1 || isempty(suffix)
    f = getfile(u,'.mesh');
else
    f = getfile(u,'.mesh',suffix);
end
