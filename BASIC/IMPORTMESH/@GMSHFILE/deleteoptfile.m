function u = deleteoptfile(u,file)
% function u = deleteoptfile(u,file)

if nargin==2
    u = setfile(u,file);
end
file = getfile(u,'.msh.opt');
if ispc
    command = 'del';
else
    command = 'rm';
end
command = [command ' ' file ';'];
dos(command);
