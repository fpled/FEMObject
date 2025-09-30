function u = compound(u,name,number)
% function u = compound(u,name,number)

if nargin<3 || isempty(number)
    number = ':';
end

u = stepcounter(u);
s = ['Compound ' name tagsintobraces(number) ';\n'];
u = addstring(u,s);
