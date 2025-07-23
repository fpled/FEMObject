function w = plus(u,v)
% function w = plus(u,v)
% addition des vecteurs u et v
%
% See also VECTOR/dot, VECTOR/norm, VECTOR/normalize, VECTOR/cross, VECTOR/abs, VECTOR/planortho, VECTOR/minus, VECTOR/rot2D,
% VECTOR/mtimes, VECTOR/times, VECTOR/mrdivide, VECTOR/ne, VECTOR/eq,
% VECTOR/plus, VECTOR/uminus, VECTOR/norm, POINT/plus, MYDOUBLEND/plus

if isa(v,'VECTOR')
    w = VECTOR(MYDOUBLEND(u)+MYDOUBLEND(v));
else
    error('')
end

