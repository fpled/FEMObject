function w = ne(u,v)
% function w = ne(u,v)
% w = u~=v
% w = 1 si la distance entre u et v est superieur a la tolerance eps
% u et v peuvent etre des multipoints
%
% See also POINT/distance, POINT/minus, POINT/mtimes, POINT/mrdivide, POINT/ne, POINT/eq,
% POINT/plus, POINT/uminus, POINT/norm, VECTOR/ne

w = (distance(u,v)>eps);    
