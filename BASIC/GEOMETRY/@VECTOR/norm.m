function w = norm(u,normtype)
% function w = norm(u,normtype)
% norme du vecteur u
% normtype : 1, 2 ou Inf (type de norme) 2 par defaut
%
% See also VECTOR/dot, VECTOR/norm, VECTOR/normalize, VECTOR/cross, VECTOR/abs, VECTOR/planortho, VECTOR/minus, VECTOR/rot2D,
% VECTOR/mtimes, VECTOR/times, VECTOR/mrdivide, VECTOR/ne, VECTOR/eq,
% VECTOR/plus, VECTOR/uminus, VECTOR/norm, POINT/norm

if nargin==1
    normtype = 2;
end

switch normtype
    case 2
        w = sqrt(sum(u.MYDOUBLEND.^2,1));
    case 1
        w = sum(abs(u.MYDOUBLEND),1);
    case Inf
        w = max(abs(u.MYDOUBLEND),[],1);
end
