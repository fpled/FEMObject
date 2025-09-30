function ctop = getctop(C)
% function ctop = getctop(C)

% Height
h = C.h;

% Top center
ctop = [0, 0, h];

% Rotation matrix
v = [C.vx, C.vy];
n = [C.nx, C.ny, C.nz];
R = calcrotation(C,v,n);

% Center
c = [C.cx, C.cy, C.cz];

% Rotate into global frame and translate to center
ctop = ctop*R + c;
