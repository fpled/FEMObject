function u = setfield(u,name,number)
% function u = setfield(u,name,number)

if nargin<3 || isempty(number)
    number = 1;
end

u = stepcounter(u);
s = ['Field[' tag2str(number) '] = ' name ';\n'];
u = addstring(u,s);
