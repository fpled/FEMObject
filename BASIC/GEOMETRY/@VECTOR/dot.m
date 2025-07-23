function w = dot(u,v)
% function w = dot(u,v)
% produit scalaire de u et v
%
% See also VECTOR/dot, VECTOR/norm, VECTOR/normalize, VECTOR/cross, VECTOR/abs, VECTOR/planortho, VECTOR/minus, VECTOR/rot2D, 
% VECTOR/mtimes, VECTOR/times, VECTOR/mrdivide, VECTOR/ne, VECTOR/eq,
% VECTOR/plus, VECTOR/uminus, VECTOR/norm

w = dot(u.MYDOUBLEND,v.MYDOUBLEND,1);
