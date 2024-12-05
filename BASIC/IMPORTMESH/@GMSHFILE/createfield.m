function u = createfield(u,name,number)
% function u = createfield(u,name,number)

if nargin<3 || isempty(number)
    number = 1;
end

u = stepcounter(u);
s = ['Field[' num2str(number) '] = ' name ' ;\n'];
u = addstring(u,s);
