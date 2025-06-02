function u = writefile(u,file)
% function u = writefile(u,file)

if nargin==2
    u = setfile(u,file);
end
file = getfile(u,'.geo');

fid = fopen(file,'w');
fprintf(fid,u.string);
fclose(fid);

u.iswritten = 1;
