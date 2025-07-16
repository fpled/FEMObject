function P = getcenter(E)
% function P = getcenter(E)

switch E.indim
    case 2
        P = [E.cx, E.cy];
    case 3
        P = [E.cx, E.cy, E.cz];
    otherwise
        error('Wrong space dimension')
end
P = POINT(P);
