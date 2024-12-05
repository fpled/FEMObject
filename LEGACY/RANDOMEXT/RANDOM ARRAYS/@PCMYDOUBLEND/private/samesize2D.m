function [u,v]=samesize2D(u,v)
utemp = u;
vtemp = v;

if isa(u,'PCMYDOUBLEND')
    u=double(u.V);
end
if isa(v,'PCMYDOUBLEND')
    v=double(v.V);
end


if size(u,1)==1 & size(u,2)==1
 u=repmat(u,[size2D(v) ones(1,ndims(u)-2)]);
elseif size(v,1)==1 & size(v,2)==1
 v=repmat(v,[size2D(u) ones(1,ndims(v)-2)]);
elseif size(u,1)==1 & (size(v,2)==size(u,2))
 u=repmat(u,[size(v,1) 1 ones(1,ndims(u)-2)]); 
elseif size(v,1)==1 & (size(v,2)==size(u,2))
 v=repmat(v,[size(u,1) 1 ones(1,ndims(u)-2)]); 
elseif size(u,2)==1 & (size(v,1)==size(u,1))
 u=repmat(u,[1 size(v,2) 1 ones(1,ndims(u)-2)]); 
elseif size(v,2)==1 & (size(v,1)==size(u,1))
 v=repmat(v,[1 size(u,2) ones(1,ndims(u)-2)]);  
elseif ~all(size2D(u)==size2D(v))
    error('les dimensions 1 � 2 doivent correspondre')
end

if isa(utemp,'PCMYDOUBLEND')
    utemp.V=MYDOUBLEND(u);
    u=utemp;
else
    u = MYDOUBLEND(u);
end
if isa(vtemp,'PCMYDOUBLEND')
    vtemp.V=MYDOUBLEND(v);
    v=vtemp;
else
    v = MYDOUBLEND(v);
end