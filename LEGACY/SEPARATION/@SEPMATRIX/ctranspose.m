function u = ctranspose(u,dim)

if nargin==1
    dim = 1:u.dim;
end

for i=1:u.m
for k=dim
    u.F{i,k}=ctranspose(u.F{i,k});
end
end

