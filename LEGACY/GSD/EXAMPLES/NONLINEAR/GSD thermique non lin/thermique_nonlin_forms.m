a = BILINFORM(1,1);
factkn=1;
g = setfact(MULTILINFORM([1,1,0,0]),factkn); % n(v,u,u,u)
n = setfact(MYFORM1(),factkn);
dn = setfact(MYFORM1GRAD(),factkn);

dim=2;
choixboun=2;
if dim==1
    nbsec =1;    
    l = LINFORM(0,10);
else
    if choixboun==1
        nbsec=1;
        l = LINFORM(0,1);   
    else
        nbsec =2;
        l = LINFORM(0,1,[],'selgroup',2);
        P = POINT([ 0,0 ; 1,0 ; 0,2 ; 1,2 ; 1,1 ; 2,1 ; 2,2 ; 0,1 ]);
        l2 = LINFORM(0,-1,[],'boundary',LINE(P(1),P(2)));
    end
end



