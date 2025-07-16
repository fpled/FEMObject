function v = gettangent(C)
% function v = gettangent(C)

v = [C.vx, C.vy];
v = v / norm(v);
v = POINT(v);
v = VECTEUR(v);
