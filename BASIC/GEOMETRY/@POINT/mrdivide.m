function u = mrdivide(u,v)
% function u = mrdivide(u,v)
% u POINT, v double
% on divise les coordonnes de u par v;
% v peut etre un MYDOUBLEND
%
% See also POINT/distance, POINT/minus, POINT/mtimes, POINT/mrdivide, POINT/ne, POINT/eq,
% POINT/plus, POINT/uminus, POINT/norm, VECTOR/mrdivide

if isa(u,'POINT')
    u.MYDOUBLEND = u.MYDOUBLEND/v;
end
