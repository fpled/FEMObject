function u = deletefile(u,ext,file)
% function u = deletefile(u,ext,file)

if nargin==3
    u = setfile(u,file);
end
if nargin<2 || isempty(ext)
    ext = '.msh.opt';
end

file = getfile(u,ext);
if ispc
    command = 'del';
else
    command = 'rm';
end
command = [command ' ' file ';'];
dos(command);
