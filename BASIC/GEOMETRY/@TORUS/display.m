function display(T)
% function display(T)

disp(' ')
disp([inputname(1) ' = (' class(T) ')'])
disp(' ')
disp(struct('cx',T.cx,'cy',T.cy,'cz',T.cz,'r1',T.r1,'r2',T.r2,'angle',T.angle))
