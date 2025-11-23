function f = getfile(u,ext,suffix)
% function f = getfile(u,ext,suffix)
% ext : character array containing file extension
% suffix : numeric array or character array containing suffix

f = u.file;
if nargin==3 && ~isempty(suffix)
    if isnumeric(suffix)
        f = [f '_' num2str(suffix)];
    else
        f = [f '_' char(suffix)];
    end
end
if nargin>=2 && ~isempty(ext) && ischar(ext)
    f = [f ext];
end
