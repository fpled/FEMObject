function lsx=eval(ls,x,cx,cy,vx,vy,funa,arga)

cx= MYDOUBLEND(reshape(full(cx),[1,1,1,numel(cx)]));
cy= MYDOUBLEND(reshape(full(cy),[1,1,1,numel(cy)]));
vx = MYDOUBLEND(reshape(full(vx),[1,1,1,numel(vx)]));
vy = MYDOUBLEND(reshape(full(vy),[1,1,1,numel(vy)]));

vn = VECTOR([vx;vy]);
vn = normalize(vn);
vt = VECTOR([vy;-vx]);
vt = normalize(vt);

delta = POINT(x)-POINT([cx,cy]);
xin = dot(delta,vn);
xit = dot(delta,vt);

if ~isa(arga,'double')
    error('l''argument de la fonction doit etre un double')
end
if size(arga,1)>1
    arga = MYDOUBLEND(reshape(full(arga'),[1,size(arga,2),1,size(arga,1)]) );
else
    arga=full(arga);
end
if isempty(arga)
    a = MYDOUBLEND(funa(double(xit)));
else
    a = MYDOUBLEND(funa(double(xit),arga));
end
lsx=xin-a;
