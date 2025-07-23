function lsx=eval(ls,x,cx,cy,a,b,vx,vy)

cx= MYDOUBLEND(reshape(full(cx),[1,1,1,numel(cx)]));
cy= MYDOUBLEND(reshape(full(cy),[1,1,1,numel(cy)]));
a = MYDOUBLEND(reshape(full(a),[1,1,1,numel(a)]));
b = MYDOUBLEND(reshape(full(b),[1,1,1,numel(b)]));
vx = MYDOUBLEND(reshape(full(vx),[1,1,1,numel(vx)]));
vy = MYDOUBLEND(reshape(full(vy),[1,1,1,numel(vy)]));

v1 = VECTOR([vx;vy]);
v1 = normalize(v1);
v2 = VECTOR([-vy;vx]);
v2 = normalize(v2);

delta = POINT(x)-POINT([cx,cy]);

xi = [dot(delta,v1)./a,dot(delta,v2)./b]; 

lsx=distance(POINT(xi),POINT([0,0]))-1;
