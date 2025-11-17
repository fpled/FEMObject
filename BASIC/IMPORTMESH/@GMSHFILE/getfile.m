function f = getfile(u,ext,iter)
% function f = getfile(u,ext,iter)
% ext : character array containing file extension
% iter : iteration number or character array containing suffix

f = u.file;
if nargin==3 && ~isempty(iter)
    if isnumeric(iter)
        f = [f '_' num2str(iter)];
    else
        f = [f '_' char(iter)];
    end
end
if nargin>=2 && ~isempty(ext) && ischar(ext)
    f = [f ext];
end
