function u = deleteoptfile(u,file)
% function u = deleteoptfile(u,file)

if nargin==2
    u = setfile(u,file);
end
file = getfilemsh(u);
file = [file '.opt'];

commande = ['rm ' file ';'];
dos(commande);