function P = getcenter(C)
% function P = getcenter(C)

switch C.indim
    case 2
        P = [C.cx, C.cy];
    case 3
        P = [C.cx, C.cy, C.cz];
    otherwise
        error('Wrong space dimension')
end
P = POINT(P);
