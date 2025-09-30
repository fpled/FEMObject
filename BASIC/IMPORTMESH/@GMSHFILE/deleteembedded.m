function u = deleteembedded(u,name,number,varargin)
% function u = deleteembedded(u,name,number)
% function u = deleteembedded(u,name,number)

if nargin<3 || isempty(number)
    number = ':';
end

u = stepcounter(u);

entities = [name tagsintobraces(number) ';'];

s = ['Delete Embedded{ ' entities ' }\n'];

u = addstring(u,s);
