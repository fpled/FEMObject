function f = getfilemsh(u,iter)
% function u = getfilemsh(u,iter)

% if nargin==1 || isempty(iter)
%     f = [u.file '.msh'];
% else
%     f = [u.file '_' num2str(iter) '.msh'];
% end

if nargin==1 || isempty(iter)
    f = getfile(u,'.msh');
else
    f = getfile(u,'.msh',iter);
end
