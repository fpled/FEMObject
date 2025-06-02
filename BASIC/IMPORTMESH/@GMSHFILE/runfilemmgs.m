function u = runfilemmgs(u,ext,options)
% function u = runfilemmgs(u,ext,options)

if nargin==1 || isempty(ext)
    ext = '.msh';
end
command = ['mmgs_O3 ' getfile(u,ext)];

if nargin==3
    command = [command ' ' options];
end

pathname = getfemobjectoptions('mmgpath');
command = fullfile(pathname,command);
dos(command);
u.ismesh = 1;
