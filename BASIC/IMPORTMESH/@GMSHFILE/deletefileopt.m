function u = deletefileopt(u,file)
% function u = deletefileopt(u,file)

% if nargin==2
%     u = setfile(u,file);
% end
% file = getfile(u,'.msh.opt');
% if ispc
%     command = 'del';
% else
%     command = 'rm';
% end
% command = [command ' ' file ';'];
% dos(command);

u = deletefile(u,'.msh.opt',file);
