function u = runfile(u,ext,options)
% function u = runfile(u,ext,options)

if nargin==1 || isempty(ext)
    ext = '.geo';
end
command = ['gmsh ' getfile(u,ext)];

if nargin<3 || isempty(options)
    options = '-format msh2';
else
    if verLessThan('matlab','9.1') % compatibility (<R2016b)
        contain = @(str,pat) ~isempty(strfind(str,pat));
    else
        contain = @contains;
    end
    if ~contain(options,'-format')
        options = [options ' -format msh2'];
    end
end
command = [command ' ' options];

pathname = getfemobjectoptions('gmshpath');
command = fullfile(pathname,command);
dos(command);
u.ismesh = 1;
