function display(C)
% function display(C)

disp(' ')
disp([inputname(1) ' = (' class(C) ')'])
disp(' ')
disp(struct('cx',C.cx,'cy',C.cy,'cz',C.cz,'r',C.r,'nx',C.nx,'ny',C.ny,'nz',C.nz,'vx',C.vx,'vy',C.vy,'h',C.h))
