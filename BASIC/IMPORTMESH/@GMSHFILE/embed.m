function u = embed(u,name,number,nameparent,numberparent)
% function u = embed(u,name,number,nameparent,numberparent)

u = stepcounter(u);
s = [name tagsintobraces(number) ' In ' nameparent '{' tag2str(numberparent) '};\n' ];
u = addstring(u,s);
