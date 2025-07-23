function w = norm(P,varargin)
% function w = norm(P,normtype)
% distance des points P a l'origine
% normtype : 1, 2 ou Inf (type de norme) 2 par defaut
%
% See also POINT/distance, POINT/minus, POINT/mtimes, POINT/mrdivide, POINT/ne, POINT/eq,
% POINT/plus, POINT/uminus, POINT/norm, VECTOR/norm

w = distance(P,0,varargin{:});
