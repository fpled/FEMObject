function var = delonlycharin(s,var)
% function var = delonlycharin(propertyname,var)
% var : tableau de cellules
% propertyname : nom d'une propriete a effacer

if isa(s,'char')
    s = {s};
end

% rep = false(1,length(s));
% pos = zeros(1,length(s));
% for j=1:length(s)
%     for i=1:length(var)
%         if isa(var{i},'char') && strcmpi(var{i},s{j})
%             rep(j) = true;
%             pos(j) = i;
%         end
%     end
% end
% 
% eff = [];
% for j=1:length(s)
%     if rep(j)
%         eff = [eff,pos(j)]; 
%     end
% end

eff = [];
for j=1:length(s)
    for i=1:length(var)
        if isa(var{i},'char') && strcmpi(var{i},s{j})
            eff = [eff, i];
        end
    end
end
eff = unique(eff); % avoid double deletion

var(eff) = [];
