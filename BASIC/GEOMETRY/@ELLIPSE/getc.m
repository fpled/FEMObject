function c = getc(E)
% function c = getc(E)

switch E.indim
    case 2
        c = [E.cx, E.cy];
    case 3
        c = [E.cx, E.cy, E.cz];
    otherwise
        error('Wrong space dimension')
end
