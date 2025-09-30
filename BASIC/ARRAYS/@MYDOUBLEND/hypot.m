function w = hypot(u,v)
% function w = hypot(u,v)

if isa(u,'MYDOUBLEND')
    w = u;
else
    w = v;
end

[u,v] = samesize(u,v);
w.double = hypot(u,v);
