function n = getnormal(C)
% function n = getnormal(C)

n = [C.nx, C.ny, C.nz];
n = n / norm(n);
n = POINT(n);
n = VECTOR(n);
