function [rep,pos] = ischarin(s,var)
% function [rep,pos] = ischarin(s,var)

if isa(s,'char')
    s = {s};
end
rep = false(1,length(s));
pos = zeros(1,length(s));
for j=1:length(s)
    for i=1:length(var)
        if isa(var{i},'char') && strcmpi(var{i},s{j})
            rep(j) = true;
            pos(j) = i;
        end
    end
end
