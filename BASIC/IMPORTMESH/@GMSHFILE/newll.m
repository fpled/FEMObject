function u = newll(u,tag)
% function u = newll(u,tag)

u = stepcounter(u);
s = [tag2str(tag) ' = newll;\n'];
u = addstring(u,s);
