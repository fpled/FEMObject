function u = createlines(u,seg,numberlines)
% function u = createlines(u,seg,numberlines)

for k=1:size(seg,1)
    if iscell(numberlines)
        u = createline(u,seg(k,:),numberlines{k});
    else
        u = createline(u,seg(k,:),numberlines(k));
    end
end
