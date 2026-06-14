function u = transfquadtri(u,number)
% function u = transfquadtri(u,number)

if nargin<2 || isempty(number)
    number = ':';
end

u = stepcounter(u);
s = ['TransfQuadTri ' tagsintobraces(number) ';\n'];
u = addstring(u,s);
