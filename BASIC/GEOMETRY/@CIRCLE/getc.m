function c = getc(C)
% function c = getc(C)

switch C.indim
    case 2
        c = [C.cx, C.cy];
    case 3
        c = [C.cx, C.cy, C.cz];
    otherwise
        error('Wrong space dimension')
end
