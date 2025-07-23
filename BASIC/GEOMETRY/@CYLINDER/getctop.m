function ctop = getctop(C)
% function ctop = getctop(C)

% Height
h = C.h;

% Top center
ctop = [0, 0, h];

%% Old version
% Rotate around axis n = [nx, ny, nz] by angle of rotation phi =
% atan2(vy, vx) using tangent vector v = [vx, vy] via Rodrigues'
% rotation formula
%% New version
% Twist the XY plane about z = [0, 0, 1] by phi = atan2(vy, vx)
% using tangent vector v = [vx, vy], then tilt from z axis to normal
% vector n = [nx, ny, nz] so that the circle's normal is n regardless
% of v = [vx, vy]
v = [C.vx, C.vy];
n = [C.nx, C.ny, C.nz];
R = calcrotation(C,v,n);

% Translate to center c = [cx, cy, cz]
c = [C.cx, C.cy, C.cz];

ctop = ctop*R + c;
