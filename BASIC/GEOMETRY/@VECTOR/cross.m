function w = cross(u,v)
% function w = cross(u,v)
% produit vectoriel de u par v
%
% See also VECTOR/dot, VECTOR/norm, VECTOR/normalize, VECTOR/cross, VECTOR/abs, VECTOR/planortho, VECTOR/minus, VECTOR/rot2D, 
% VECTOR/mtimes, VECTOR/times, VECTOR/mrdivide, VECTOR/ne, VECTOR/eq,
% VECTOR/plus, VECTOR/uminus, VECTOR/norm, MYDOUBLEND/cross

w = VECTOR(cross(MYDOUBLEND(u),MYDOUBLEND(v),1));
