function u = news(u,tag)
% function u = news(u,tag)

u = stepcounter(u);
s = [tag2str(tag) ' = news;\n'];
u = addstring(u,s);
