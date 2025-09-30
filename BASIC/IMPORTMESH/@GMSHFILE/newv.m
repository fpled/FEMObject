function u = newv(u,tag)
% function u = newv(u,tag)

u = stepcounter(u);
s = [tag2str(tag) ' = newv;\n'];
u = addstring(u,s);
