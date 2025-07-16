function u = hypot(u,v)
% function u = hypot(u,v)

if isa(u,'MYDOUBLEND')
    w = u;
else
    w = v;
end

[u,v] = samesize(u,v);
w.double = hypot(u,v);
