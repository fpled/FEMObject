function u = newcl(u,tag)
% function u = newcl(u,tag)

u = stepcounter(u);
s = [tag2str(tag) ' = newcl;\n'];
u = addstring(u,s);
