function u = rot2D(u,a)
% function u = rot2D(u,a)
% rotation du vecteur u d'un angle a
% marche uniquement en 2D
%
% See also VECTOR/dot, VECTOR/norm, VECTOR/normalize, VECTOR/cross, VECTOR/abs, VECTOR/planortho, VECTOR/minus, VECTOR/rot2D, 
% VECTOR/mtimes, VECTOR/times, VECTOR/mrdivide, VECTOR/ne, VECTOR/eq,
% VECTOR/plus, VECTOR/uminus, VECTOR/norm

R = [cos(a),-sin(a);sin(a),cos(a)];

u.MYDOUBLEND = R*u.MYDOUBLEND;
