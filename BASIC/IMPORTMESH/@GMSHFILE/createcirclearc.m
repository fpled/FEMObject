function u = createcirclearc(u,numcenter,numpoints,numbercurve)
% function u = createcirclearc(u,numcenter,numpoints,numbercurve)

if length(numpoints)~=2
    error('A circle arc is defined by the 2 indices of the start and end points: numpoints = [start,end]')
end

u = createentity(u,'Circle',[numpoints(1),numcenter,numpoints(2)],numbercurve);
