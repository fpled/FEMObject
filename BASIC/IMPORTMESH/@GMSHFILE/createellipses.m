function u = createellipses(u,center,radii,angles,numbercurves)
% function u = createellipses(u,center,radii,angles,numbercurves)

for k=1:size(angles,1)
    u = createellipse(u,center,radii,numbercurves(k),angles(k,:));
end
