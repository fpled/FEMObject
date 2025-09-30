function u = setmaxtag(u,name,number)
% function u = setmaxtag(u,name,number)

if nargin<2 || isempty(name)
    name = 'Point';
end
if nargin<3 || isempty(number)
    number = 1;
end
if ~ischar(number)
    number = num2str(number);
end

u = stepcounter(u);
s = ['SetMaxTag ' name '(' number ');\n'];
u = addstring(u,s);
