function f = getfile(u,ext,iter)
% function f = getfile(u,ext,iter)
% ext : character array containing file extension
% iter : iteration number

f = u.file;
if nargin==3 && ~isempty(iter) && isnumeric(iter)
    f = [f '_' num2str(iter)];
end
if nargin>=2 && ~isempty(ext) && ischar(ext)
    f = [f ext];
end
