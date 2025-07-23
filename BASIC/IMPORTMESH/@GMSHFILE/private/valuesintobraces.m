function s = valuesintobraces(values)
% function s = valuesintobraces(values)

% s = [' {' num2str(values(1))];
% for k=2:length(values)
%     s = [ s, ' , '  num2str(values(k)) ];
% end
% s = [s, '} '];

if isempty(values)
    s = ' {} ';
    return;
end

if ischar(values) && isrow(values)
    vals = {values};
elseif ischar(values) && ismatrix(values)
    vals = cellstr(values);
elseif isstring(values)
    vals = cellstr(values);
elseif isnumeric(values)
    vals = arrayfun(@num2str,values(:)','UniformOutput',false);
elseif iscell(values)
    vals = cell(size(values(:)'));
    for k=1:numel(values)
        v = values{k};
        if isnumeric(v)
            vals{k} = num2str(v);
        elseif isstring(v)
            vals{k} = char(v);
        elseif ischar(v)
            vals{k} = v;
        else
            error('Unsupported element type in cell array.');
        end
    end
else
    error('Unsupported input type for valuesintobraces.');
end

nums = strjoin(vals, ' , ');
s = [' {' nums '} '];

