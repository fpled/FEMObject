function s = tagsintobraces(tags)
% function s = tagsintobraces(tags)
% Convert mixed input tags into a list of char tags within braces for Gmsh.

% s = ['{' num2str(tags(1))];
% for k=2:length(tags)
%     s = [ s, ', '  num2str(tags(k)) ];
% end
% s = [s, '}'];

if isempty(tags)
    s = '{}';
    return;
end

if iscell(tags)
    T = cellfun(@tag2str, tags(:).', 'UniformOutput', false);
    
elseif isnumeric(tags) || islogical(tags)
    t = tags(:).';
    if ~isreal(t)
        error('tagsintobraces:NumericMustBeReal', ...
              'Numeric tags must be real (no complex).');
    end
    if ~all(isfinite(t))
        error('tagsintobraces:NumericMustBeFinite', ...
              'Numeric tags must be finite (no NaN/Inf).');
    end
    T = arrayfun(@num2str, t, 'UniformOutput', false);
    
elseif isstring(tags)
    t = tags(:).';
    if any(ismissing(t))
        error('tagsintobraces:StringMissingNotAllowed', ...
              'Missing strings are not allowed in tag lists.');
    end
    T = arrayfun(@tag2str, t, 'UniformOutput', false);
    
elseif ischar(tags)
    if isrow(tags)
        T = {tags};
    else
        T = cellstr(tags).';
    end
    
else
    error('tagsintobraces:UnsupportedType', ...
          'Unsupported input type "%s".', class(tags));
end

s = ['{' strjoin(T, ', ') '}'];

end
