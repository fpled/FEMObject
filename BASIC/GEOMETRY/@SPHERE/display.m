function display(S)
% function display(S)

disp(' ')
disp([inputname(1) ' = (' class(S) ')'])
disp(' ')
disp(struct('cx',S.cx,'cy',S.cy,'cz',S.cz,'r',S.r,'nx',S.nx,'ny',S.ny,'nz',S.nz,'vx',S.vx,'vy',S.vy))
