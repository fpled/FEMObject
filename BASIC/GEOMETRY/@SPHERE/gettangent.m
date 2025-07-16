function v = gettangent(S)
% function v = gettangent(S)

v = [S.vx, S.vy];
v = v / norm(v);
v = POINT(v);
v = VECTEUR(v);
