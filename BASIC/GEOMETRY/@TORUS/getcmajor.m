function cmajor = getcmajor(T)
% function cmajor = getcmajor(T)

% Minor radius
r2 = T.r2;

cmajor = [0, 0,  r2;  % +z
          0, 0, -r2]; % -z

% Rotation matrix
v = [T.vx, T.vy];
n = [T.nx, T.ny, T.nz];
R = calcrotation(T,v,n);

% Center
c = [T.cx, T.cy, T.cz];

% Rotate into global frame and translate to center
cmajor = cmajor*R + c;

end
