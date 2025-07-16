function n = getnormal(S)
% function n = getnormal(S)

n = [S.nx, S.ny, S.nz];
n = n / norm(n);
n = POINT(n);
n = VECTEUR(n);
