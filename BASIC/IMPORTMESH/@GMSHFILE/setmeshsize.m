function u = setmeshsize(u,value,number)
% function u = setmeshsize(u,value,number)

if nargin<3 || isempty(number)
    number = ':';
end

u = stepcounter(u);
s = ['MeshSize ' valuesintobraces(number) ' = ' num2str(value) ' ;\n'];
u = addstring(u,s);
