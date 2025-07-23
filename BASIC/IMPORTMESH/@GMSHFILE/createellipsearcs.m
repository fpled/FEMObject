function u = createellipsearcs(u,numbercenter,seg,numbermajorpoints,numbercurves)
% function u = createellipsearcs(u,numbercenter,seg,numbermajorpoints,numbercurves)

for k=1:size(seg,1)
    u = createellipsearc(u,numbercenter,seg(k,:),numbermajorpoints(k),numbercurves(k));
end
