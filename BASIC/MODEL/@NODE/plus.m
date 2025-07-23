function u = plus(u,v)

if isa(u,'NODE') && isa(v,'VECTOR')
    u.POINT = u.POINT + v ; 
elseif isa(v,'NODE') && isa(u,'VECTOR')
    v.POINT = v.POINT + u ; 
end

    