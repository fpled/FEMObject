function u = openfile(u,ext)
% function u = openfile(u,ext)

if nargin==1 || isempty(ext)
    ext = '.geo';
end
command = ['gmsh ' getfile(u,ext)];

pathname = getfemobjectoptions('gmshpath');
command = fullfile(pathname,command);
dos(command);
