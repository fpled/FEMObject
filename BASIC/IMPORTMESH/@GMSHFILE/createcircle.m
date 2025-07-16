function u = createcircle(u,numcenter,numpoints,numberline)
% function u = createcircle(u,numcenter,numpoints,numberline)

if length(numpoints)~=2
    error('pour creer un cercle, il faut 2 numeros de points')
end
u = createentity(u,'Circle',[numpoints(1),numcenter,numpoints(2)],numberline);
