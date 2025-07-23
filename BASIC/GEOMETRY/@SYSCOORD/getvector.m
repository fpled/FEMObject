function v = getvector(v,dim)
% function v = getvector(v,dim)

if nargin==1
    v = VECTOR(v.MYDOUBLEND);
else
    v = VECTOR(v.MYDOUBLEND(:,dim));
end
