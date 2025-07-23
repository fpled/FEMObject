function w = times(u,v)
% function w = mtimes(u,v)
% on multiplie composante par composante u et v
% on peut egalement entrer des MYDOUBLEND
%
% See also VECTOR/dot, VECTOR/norm, VECTOR/normalize, VECTOR/cross, VECTOR/abs, VECTOR/planortho, VECTOR/minus, VECTOR/rot2D, 
% VECTOR/mtimes, VECTOR/times, VECTOR/mrdivide, VECTOR/ne, VECTOR/eq,
% VECTOR/plus, VECTOR/uminus, VECTOR/norm, MYDOUBLEND/times

w = dot(MYDOUBLEND(u),MYDOUBLEND(v),1);
