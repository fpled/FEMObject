function n = getnormal(T)
% function n = getnormal(T)

n = [T.nx, T.ny, T.nz];
n = n / norm(n);
n = POINT(n);
n = VECTOR(n);
