function u = refinemesh(u)
% function u = refinemesh(u)

u = stepcounter(u);
s = 'RefineMesh;\n';
u = addstring(u,s);
