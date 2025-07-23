function w = eq(u,v)
% function w = eq(u,v)
% w = u==v
% w = 1 si la distance entre u et v est inferieure a la tolerance eps
% u et v peuvent etre des multipoints
%
% See also POINT/distance, POINT/minus, POINT/mtimes, POINT/mrdivide, POINT/ne, POINT/eq,
% POINT/plus, POINT/uminus, POINT/norm, VECTOR/eq

w = (distance(u,v)<=eps);    
