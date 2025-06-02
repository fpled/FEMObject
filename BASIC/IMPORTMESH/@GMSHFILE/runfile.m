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
        noFormat = isempty(strfind(options,'-format'));
    else
        noFormat = ~contains(options,'-format');
    end
    if noFormat
        options = [options ' -format msh2'];
    end
end
command = [command ' ' options];

pathname = getfemobjectoptions('gmshpath');
command = fullfile(pathname,command);
dos(command);
u.ismesh = 1;
