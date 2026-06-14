function u = recombinemesh(u)
% function u = recombinemesh(u)

u = stepcounter(u);
s = 'RecombineMesh;\n';
u = addstring(u,s);
