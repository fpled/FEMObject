function f = getfilesol(u,iter)
% function u = getfilesol(u,iter)

if nargin==1 || isempty(iter)
    f = [u.file '.sol'];
else
    f = [u.file '_' num2str(iter) '.sol'];
end
