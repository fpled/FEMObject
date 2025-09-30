function u = recombinesurface(u,number)
% function u = recombinesurface(u,number)

if nargin<2 || isempty(number)
    number = ':';
end

u = stepcounter(u);
s = ['Recombine Surface' tagsintobraces(number) ';\n'];
u = addstring(u,s);
