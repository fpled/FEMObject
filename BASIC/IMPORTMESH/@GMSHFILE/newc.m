function u = newc(u,tag)
% function u = newc(u,tag)

u = stepcounter(u);
s = [tag2str(tag) ' = newc;\n'];
u = addstring(u,s);
