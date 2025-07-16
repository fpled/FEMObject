function v = gettangent(E)
% function v = gettangent(E)

v = [E.vx, E.vy];
v = v / norm(v);
v = POINT(v);
v = VECTEUR(v);
