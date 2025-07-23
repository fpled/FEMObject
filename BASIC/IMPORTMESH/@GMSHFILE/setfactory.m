function u = setfactory(u,name)
% function u = setfactory(u,name)

if nargin<2 || isempty(name)
    name = "Built-in";
end

u = stepcounter(u);
s = ['SetFactory("' char(name) '");\n'];
u = addstring(u,s);
