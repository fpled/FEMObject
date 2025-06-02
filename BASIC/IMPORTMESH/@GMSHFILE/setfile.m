function u = setfile(u,file)
% function u = setfile(u,file)

[filepath,name,ext] = fileparts(file);

u.file = fullfile(filepath,name);

if strcmpi(ext,'.geo')
    u.iswritten = 1;
elseif strcmpi(ext,'.msh') || strcmpi(ext,'.mesh')
    u.ismesh = 1;
end
