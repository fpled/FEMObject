function u = mrdivide(u,v)
% function u = mrdivide(u,v)
% on divise les composantes de u par v (double ou MYDOUBLEND)
%
% See also VECTOR/dot, VECTOR/norm, VECTOR/normalize, VECTOR/cross, VECTOR/abs, VECTOR/planortho, VECTOR/minus, VECTOR/rot2D,
% VECTOR/mtimes, VECTOR/times, VECTOR/mrdivide, VECTOR/ne, VECTOR/eq,
% VECTOR/plus, VECTOR/uminus, VECTOR/norm, MYDOUBLEND/mrdivide

if (isa(u,'VECTOR') )
    u.MYDOUBLEND = u.MYDOUBLEND/v;
else
    error('operation non definie')
end
