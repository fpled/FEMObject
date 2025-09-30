function v = gettangent(T)
% function v = gettangent(T)

v = [T.vx, T.vy];
v = v / norm(v);
v = POINT(v);
v = VECTOR(v);
