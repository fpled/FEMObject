function w = minus(u,v)
% function w = minus(u,v)
% si u POINT et v POINT : w = vu (VECTOR)
% sinon w = VECTOR(u)-VECTOR(v)
%
% See also POINT/distance, POINT/minus, POINT/mtimes, POINT/mrdivide, POINT/ne, POINT/eq,
% POINT/plus, POINT/uminus, POINT/norm, VECTOR/minus


if isa(u,'POINT') && isa(v,'POINT')
    w = VECTOR(u)-VECTOR(v);
else
    w = plus(u,uminus(v));
end