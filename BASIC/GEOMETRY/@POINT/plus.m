function w = plus(u,v)
% function w = plus(u,v)
% si u POINT et v POINT : w = Ou + Ov
% si u POINT et v VECTEUR : w = Ou+v;
% sinon si u POINT : w = Ou +v
%       si v POINT : w = u +Ov
%
% See also POINT/distance, POINT/minus, POINT/mtimes, POINT/mrdivide, POINT/ne, POINT/eq,
% POINT/plus, POINT/uminus, POINT/norm, VECTEUR/plus

if isa(v,'VECTEUR')
    w = POINT(VECTEUR(u)+v);
elseif isa(u,'POINT')
    w = POINT(MYDOUBLEND(u)+double(v));
elseif isa(v,'POINT')
    w = POINT(MYDOUBLEND(v)+double(u));
end