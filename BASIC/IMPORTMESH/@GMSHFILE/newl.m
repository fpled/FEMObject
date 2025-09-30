function u = newl(u,tag)
% function u = newl(u,tag)

u = stepcounter(u);
s = [tag2str(tag) ' = newl;\n'];
u = addstring(u,s);
