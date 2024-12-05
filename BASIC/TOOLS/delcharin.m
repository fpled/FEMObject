function var = delcharin(s,var,def)
% function var = delcharin(propertyname,var)
% var : tableau de cellules contenant des paires du type {'propertyname',propertyvalue}
% propertyname : nom d'une propriete a effacer

if isa(s,'char')
    s={s};
end

rep = zeros(1,length(s));
pos = zeros(1,length(s));
for j=1:length(s)
    for i=1:length(var)
        if isa(var{i},'char') && strcmpi(var{i},s{j})
            rep(j)=1;
            pos(j)=i;
        end
    end
end

eff = [];
for j=1:length(s)
    if rep(j)  
        eff = [eff,pos(j),pos(j)+1]; 
    end
end

var(eff)=[];
