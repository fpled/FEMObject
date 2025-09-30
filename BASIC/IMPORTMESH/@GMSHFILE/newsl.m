function u = newsl(u,tag)
% function u = newsl(u,tag)

u = stepcounter(u);
s = [tag2str(tag) ' = newsl;\n'];
u = addstring(u,s);
