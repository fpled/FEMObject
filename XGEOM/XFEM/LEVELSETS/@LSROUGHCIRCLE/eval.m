function lsx=eval(ls,x,cx,cy,vx,vy,funr,argr)

cx= MYDOUBLEND(reshape(full(cx),[1,1,1,numel(cx)]));
cy= MYDOUBLEND(reshape(full(cy),[1,1,1,numel(cy)]));
vx= MYDOUBLEND(reshape(full(vx),[1,1,1,numel(vx)]));
vy= MYDOUBLEND(reshape(full(vy),[1,1,1,numel(vy)]));

v = VECTOR([vx;vy]);
v2 = rot2D(v,pi/2);

c = POINT([cx,cy]);
X = dot(POINT(x)-c,v);
Y = dot(POINT(x)-c,v2);
theta = atan2(double(Y),double(X));

if ~isa(argr,'double')
    error('l''argument de la fonction doit etre un double')
end
if size(argr,1)>1
    argr = MYDOUBLEND(reshape(full(argr'),[1,size(argr,2),1,size(argr,1)]) );
end
if isempty(argr)
    r = funr(theta);
else
    r = funr(theta,argr);
end
lsx=distance(POINT(x),c)-r;
