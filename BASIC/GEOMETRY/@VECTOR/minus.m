function u = minus(u,v)
% function u = minus(u,v)
% soustraction des vecteurs u et v
%
% See also VECTOR/dot, VECTOR/norm, VECTOR/normalize, VECTOR/cross, VECTOR/abs, VECTOR/planortho, VECTOR/minus, VECTOR/rot2D, 
% VECTOR/mtimes, VECTOR/times, VECTOR/mrdivide, VECTOR/ne, VECTOR/eq,
% VECTOR/plus, VECTOR/uminus, VECTOR/norm, POINT/minus, MYDOUBLEND/minus

u.MYDOUBLEND =  u.MYDOUBLEND - v.MYDOUBLEND ;
