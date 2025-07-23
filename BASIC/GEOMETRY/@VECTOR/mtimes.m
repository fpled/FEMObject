function w = mtimes(u,v)
% function w = mtimes(u,v)
% on multiplie u par v
% si v VECTOR et u double (ou MYDOUBLEND) 
% alors u doit verifier size(u,2)=1 ou size(u,2)=getindim(v) et size(u,2)=getindim(v) 
% si u VECTOR et v double alors v doit verifier size(v,1)=1
% on peut egalement entrer des MYDOUBLEND
%
% See also VECTOR/dot, VECTOR/norm, VECTOR/normalize, VECTOR/cross, VECTOR/abs, VECTOR/planortho, VECTOR/minus, VECTOR/rot2D, 
% VECTOR/mtimes, VECTOR/times, VECTOR/mrdivide, VECTOR/ne, VECTOR/eq,
% VECTOR/plus, VECTOR/uminus, VECTOR/norm, MYDOUBLEND/mtimes

w = VECTOR(MYDOUBLEND(u)*MYDOUBLEND(v));
