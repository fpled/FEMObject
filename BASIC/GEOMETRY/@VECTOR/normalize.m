function u = normalize(u,varargin)
% function u = normalize(u,normtype)
% normalisation du vecteur u 
% normtype : 1, 2 ou Inf (type de norme) 2 par defaut
%
% See also VECTOR/dot, VECTOR/norm, VECTOR/normalize, VECTOR/cross, VECTOR/abs, VECTOR/planortho, VECTOR/minus, VECTOR/rot2D, 
% VECTOR/mtimes, VECTOR/times, VECTOR/mrdivide, VECTOR/ne, VECTOR/eq,
% VECTOR/plus, VECTOR/uminus, VECTOR/norm

n = norm(u,varargin{:});
u.MYDOUBLEND = u.MYDOUBLEND./n;
