function u = setcoordinatessurface(u,name)
% function u = setcoordinatessurface(u,name)

u = stepcounter(u);
s = ['Coordinates Surface ' name ';\n'];
u = addstring(u,s);
