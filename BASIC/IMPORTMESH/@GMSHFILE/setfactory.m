function u = setfactory(u,name)
% function u = setfactory(u,name)

if nargin<2 || isempty(name)
    name = 'Built-in';
end
if isstring(name)
    name = char(name);
end

u = stepcounter(u);
s = ['SetFactory("' name '");\n'];
u = addstring(u,s);
