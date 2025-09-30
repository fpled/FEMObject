function u = newreg(u,tag)
% function u = newreg(u,tag)

u = stepcounter(u);
s = [tag2str(tag) ' = newreg;\n'];
u = addstring(u,s);
