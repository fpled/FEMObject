function u = createcircles(u,center,radius,angles,numbercurves)
% function u = createcircles(u,center,radius,angles,numbercurves)

for k=1:size(angles,1)
    u = createcircle(u,center,radius,numbercurves(k),angles(k,:));
end
