function f = getfilemesh(u,iter)
% function u = getfilemesh(u,iter)

% if nargin==1 || isempty(iter)
%     f = [u.file '.mesh'];
% else
%     f = [u.file '_' num2str(iter) '.mesh'];
% end

if nargin==1 || isempty(iter)
    f = getfile(u,'.mesh');
else
    f = getfile(u,'.mesh',iter);
end
