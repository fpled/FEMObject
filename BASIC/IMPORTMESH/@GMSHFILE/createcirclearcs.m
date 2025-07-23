function u = createcirclearcs(u,numbercenter,seg,numbercurves)
% function u = createcirclearcs(u,numbercenter,seg,numbercurves)

for k=1:size(seg,1)
    u = createcirclearc(u,numbercenter,seg(k,:),numbercurves(k));
end
