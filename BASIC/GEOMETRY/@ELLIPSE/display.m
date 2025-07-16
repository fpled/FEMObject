function display(E)
% function display(E)

disp(' ')
disp([inputname(1) ' = (' class(E) ')'])
disp(' ')
switch E.indim
    case 2
        disp(struct('cx',E.cx,'cy',E.cy,'a',E.a,'b',E.b,'vx',E.vx,'vy',E.vy))
    case 3
        disp(struct('cx',E.cx,'cy',E.cy,'cz',E.cz,'a',E.a,'b',E.b,'nx',E.nx,'ny',E.ny,'nz',E.nz,'vx',E.vx,'vy',E.vy))
end
