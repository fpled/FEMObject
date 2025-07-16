function display(C)
% function display(C)

disp(' ')
disp([inputname(1) ' = (' class(C) ')'])
disp(' ')
switch C.indim
    case 2
        disp(struct('cx',C.cx,'cy',C.cy,'a',C.a,'b',C.b,'vx',C.vx,'vy',C.vy))
    case 3
        disp(struct('cx',C.cx,'cy',C.cy,'cz',C.cz,'a',C.a,'b',C.b,'nx',C.nx,'ny',C.ny,'nz',C.nz,'vx',C.vx,'vy',C.vy))
end
