function w = eq(u,v)
% function w = eq(u,v)
% w = u==v
% w = 1 si les vecteurs sont confondus a la tolerance eps pres
% u et v peuvent etre des multivecteurs
%
% See also VECTOR/dot, VECTOR/norm, VECTOR/normalize, VECTOR/cross, VECTOR/abs, VECTOR/planortho, VECTOR/minus, VECTOR/rot2D, 
% VECTOR/mtimes, VECTOR/times, VECTOR/mrdivide, VECTOR/ne, VECTOR/eq,
% VECTOR/plus, VECTOR/uminus, VECTOR/norm, POINT/eq

w = (norm(u-v)<eps);
