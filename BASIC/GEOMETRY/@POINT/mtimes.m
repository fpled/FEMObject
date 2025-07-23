function w = mtimes(u,v)
% function w = mtimes(u,v)
% si u POINT, v double (ou MYDOUBLEND)
% on multiplie les coordonnes de u par v
% si v POINT, u double (ou MYDOUBLEND)
% on multiplie les coordonnes de v par u
%
% See also POINT/distance, POINT/minus, POINT/mtimes, POINT/mrdivide, POINT/ne, POINT/eq,
% POINT/plus, POINT/uminus, POINT/norm, VECTOR/mtimes

w = POINT(MYDOUBLEND(u)*MYDOUBLEND(v));
