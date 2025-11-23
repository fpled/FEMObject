function f = getfilemsh(u,suffix)
% function u = getfilemsh(u,suffix)

if nargin==1 || isempty(suffix)
    f = getfile(u,'.msh');
else
    f = getfile(u,'.msh',suffix);
end
