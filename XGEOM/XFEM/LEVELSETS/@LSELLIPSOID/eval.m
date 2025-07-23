function lsx=eval(ls,x,cx,cy,cz,a,b,c,vax,vay,vaz,vbx,vby,vbz)

cx= MYDOUBLEND(reshape(full(cx),[1,1,1,numel(cx)]));
cy= MYDOUBLEND(reshape(full(cy),[1,1,1,numel(cy)]));
cz= MYDOUBLEND(reshape(full(cz),[1,1,1,numel(cz)]));
a = MYDOUBLEND(reshape(full(a),[1,1,1,numel(a)]));
b = MYDOUBLEND(reshape(full(b),[1,1,1,numel(b)]));
c = MYDOUBLEND(reshape(full(c),[1,1,1,numel(c)]));
vax = MYDOUBLEND(reshape(full(vax),[1,1,1,numel(vax)]));
vay = MYDOUBLEND(reshape(full(vay),[1,1,1,numel(vay)]));
vaz = MYDOUBLEND(reshape(full(vaz),[1,1,1,numel(vaz)]));
vbx = MYDOUBLEND(reshape(full(vbx),[1,1,1,numel(vbx)]));
vby = MYDOUBLEND(reshape(full(vby),[1,1,1,numel(vby)]));
vbz = MYDOUBLEND(reshape(full(vbz),[1,1,1,numel(vbz)]));

va = VECTOR([vax;vay;vaz]);
va = normalize(va);
vb = VECTOR([vbx;vby;vbz]);
vb = normalize(vb);
vc = cross(va,vb);

delta = POINT(x)-POINT([cx,cy,cz]);

xi = [dot(delta,va)./a,dot(delta,vb)./b,dot(delta,vc)./c]; 

lsx=distance(POINT(xi),POINT([0,0,0]))-1;
