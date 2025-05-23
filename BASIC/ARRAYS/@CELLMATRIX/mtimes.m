function w=mtimes(u,v,masse,n)

if issparse(u) & ~issparse(v)
    v=sparse(v);
elseif issparse(v) & ~issparse(u)
    u=sparse(u);  
end
    

if isa(u,'MULTIMATRIX') & isa(v,'MULTIMATRIX') & all(u.sm==v.sm)
if (length(u.s)>2 && any(u.s(3:end)>1)) || (length(v.s)>2 && any(v.s(3:end)>1))
if length(u.s)==2 || length(v.s)==2
    error('les deux MULTIMATRIX soivent avoir la meme 3eme dimension')
end
su2 = u.s(1:2);
sv2 = v.s(1:2);
su3 = u.s(3:end);
sv3 = v.s(3:end);

error('pas programme');

else
    if all(u.s==1) & all(v.s==1)
        
    if nargin==2
        w = u ;
        w.value = times(u.value,v.value);
        
    elseif (nargin==3 | (nargin==4 && n==1)) & isa(masse,'MULTIMATRIX')
        
        w = u ;
        vtemp = double(v.value*masse);
        w.value = u.value*vtemp';
        
    elseif nargin==4 & isa(masse,'MULTIMATRIX') & n==2
        w = v ;
        vtemp = double(u.value*masse);
        w.value = v.value*vtemp';
        
    else
        error('non programme')
        
    end
       
        
    elseif all(u.s==1) | all(v.s==1)

        
    if all(u.s==1)
     u=repmat(u,v.s);
    else
     v=repmat(v,u.s);   
    end
    
    if nargin==2
    w=times(u,v);  
    elseif nargin==3
    w=times(u,v,masse);      
    elseif nargin==4
    w=times(u,v,masse,n);  
    end
    
    
    else
    s = [u.s(1),v.s(2)];
    value = [] ;
    i=[];
    j=[];
    
    if nargin==2        
    p=size(u.value,2);
    
    for k=1:p
    vtemp = reshape(v.value(:,k),v.s);
    utemp = reshape(u.value(:,k),u.s);
    uvtemp = utemp * vtemp ; 
    Ik=find(uvtemp);Ik=Ik(:);
    Jk=k*ones(length(Ik),1);    
    i=[i;Ik];
    j=[j;Jk];
    uvtemp = uvtemp(Ik);
    value=[value;uvtemp(:)];
    end
    elseif isa(masse,'MULTIMATRIX') & (nargin==3 | n==1)
    p=size(v.value,2);
    multivtemp = v.value*masse ;
    for k=1:size(u.value,2)
    vtemp = getmatrix(multivtemp,k);
    vtemp = reshape(vtemp,[v.s(1), v.s(2)*size(vtemp,2)]);
    utemp = reshape(u.value(:,k),u.s);
    uvtemp = utemp * vtemp ; 
    uvtemp =reshape(uvtemp,[u.s(1)*v.s(2),size(v.value,2)]);
    [ik,jk]=find(uvtemp);
    i=[i;ik(:)];
    j=[j;jk(:)];
    Ik = find(uvtemp);
    uvtemp = uvtemp(Ik);
    value=[value;uvtemp(:)];
    end
    elseif isa(masse,'MULTIMATRIX') & n==2
    p=size(u.value,2);
    
    multiutemp = u.value*masse ;
    for k=1:size(v.value,2)    
    utemp = getmatrix(multiutemp,k);
    utemp = reshape(utemp.',[size(utemp,2)*u.s(1),u.s(2)]);
    vtemp = reshape(v.value(:,k),v.s);    
    uvtemp = utemp * vtemp ; 

    uvtemp =reshape(uvtemp,[size(u.value,2),u.s(1)*v.s(2)]).';
    [ik,jk]=find(uvtemp);
    i=[i;ik(:)];
    j=[j;jk(:)];
    Ik = find(uvtemp);
    uvtemp = uvtemp(Ik);
    value=[value;uvtemp(:)];
    end    
    end
    
    value = sparse(i,j,value,prod(s),p);
    
    w = u ; 
    w.value = value ; 
    w.s = s;
    end
    
end

elseif isa(u,'MULTIMATRIX') & isa(v,'double')
    
    if all(size(v)==1)
    s=u.s;
    value = u.value*v;
    else    
    p=size(u.value,2);
    s = [u.s(1),size(v,2)];
    utemp  = reshape(u.value.',p*u.s(1),u.s(2));
    value = reshape(utemp * v,p,u.s(1)*size(v,2)).' ;
    end
    
    w = u ; 
    w.value = value ; 
    w.s = s;
    
elseif isa(v,'MULTIMATRIX') & isa(u,'double')

    if all(size(u)==1)
    s=v.s;
    value = u*v.value;
    else    
    p=size(v.value,2);
    s = [size(u,1),v.s(2)];
    vtemp  = reshape(v.value,v.s(1),v.s(2)*p);
    value = reshape(u * vtemp,size(u,1)*v.s(2),p) ;
    end

    w = v ; 
    w.value = value ; 
    w.s = s;
    
elseif isa(u,'MULTIMATRIX') & isa(v,'MULTIMATRIX') & ~all(u.sm==v.sm)
if length(u.sm)>2 | length(v.sm)>2
    error('pas programme')
end
if all(u.sm==1)
    w = reshape(u.value,u.s)*v
elseif all(v.sm==1)

w = u*reshape(v.value,v.s);

elseif u.sm(1)==1 & v.sm(2)==1
v.value = repmat(v.value,1,u.sm(2));
u.value = repmat(u.value,1,v.sm(1));
sm = [v.sm(1),u.sm(2)];
u.sm=sm;
v.sm=sm;
u = multitranspose(u);
w = mtimes(u,v);
    
elseif u.sm(2)==1 & v.sm(1)==1

v.value = repmat(v.value,1,u.sm(1));
u.value = repmat(u.value,1,v.sm(2));
sm = [u.sm(1),v.sm(2)];
u.sm=sm;
v.sm=sm;
v=multitranspose(v);
w = mtimes(u,v);

elseif u.sm(1)==1 & u.sm(2)==v.sm(2)
u.value = repmat(u.value,1,v.sm(1));
sm = [v.sm(1),u.sm(2)];
u.sm=sm;
u = multitranspose(u);
w = mtimes(u,v);

elseif u.sm(2)==1 & v.sm(1)==u.sm(1)

u.value = repmat(u.value,1,v.sm(2));
sm = [u.sm(1),v.sm(2)];
u.sm=sm;
w = mtimes(u,v);

elseif u.sm(1)==v.sm(1) & v.sm(2)==1
v.value = repmat(v.value,1,u.sm(2));
sm = [v.sm(1),u.sm(2)];
v.sm=sm;
w = mtimes(u,v);
    
elseif u.sm(2)==v.sm(2) & v.sm(1)==1

v.value = repmat(v.value,1,u.sm(1));
sm = [u.sm(1),v.sm(2)];
v.sm=sm;
v=multitranspose(v);
w = mtimes(u,v);
else
 
    error('pas defini')
end


else
       error('mtimes not defined') 
end

