function u = newp(u,tag)
% function u = newp(u,tag)

u = stepcounter(u);
s = [tag2str(tag) ' = newp;\n'];
u = addstring(u,s);
