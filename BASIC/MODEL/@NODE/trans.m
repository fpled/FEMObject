function N = trans(N,v,liste)

if isa(v,'VECTOR')
    N.POINT=N.POINT+v;
end