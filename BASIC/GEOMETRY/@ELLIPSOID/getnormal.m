function n = getnormal(E)
% function n = getnormal(E)

n = [E.nx, E.ny, E.nz];
n = n / norm(n);
n = POINT(n);
n = VECTEUR(n);
