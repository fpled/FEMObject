%% modele stochastique
r=10;
p=8;

P = POINT([ 0,0 ; 1,0 ; 0,2 ; 1,2 ; 1,1 ; 2,1 ; 2,2 ]);
S1 = mesh(DOMAIN(2,P(1),P(4)),r,r*2);
S2 = mesh(DOMAIN(2,P(5),P(7)),r,r);
S = union(S1,S2);
S=concatgroupelem(S);
S=convertelem(S,'TRI3');
secdet=0;

%RV = RVUNIFORM(1/1000,1/100);
RV = RVUNIFORM(1/10000,1/100)
X = PCMODEL(RV,'pcg','order',p);
PC = getPC(X);

k=X{1};
cx=0;
cy=1;


mat=FOUR_ISOT('k',k,'b',[cx;cy]);
S = setmaterial(S,mat);
S=final(S);
L1 = LINE(P(6),P(7));
L2 = LINE(P(1),P(3));
L3 = LINE(P(3),P(7));
L4 = LINE(P(6),P(5));
L5 = LINE(P(2),P(5));
L6 = LINE(P(1),P(2));

% conditions aux limites
S=addcl(S,L1,'T',0);
S=addcl(S,L2,'T',0);
S=addcl(S,L3,'T',0);
S=addcl(S,L4,'T',0);
S=addcl(S,L5,'T',0);
S=addcl(S,L6,'T',0);

f=bodyloadPG(S,[],'QN',1);

Kpc = calc_rigipc(S,PC);   
fpc = f*one(PC);
upc = solve(Kpc,fpc);

nbnode = S.nbnode;

mat=FOUR_ISOT('k',1,'b',[cx;cy]);
S = setmaterial(S,mat);
Kdiff1 = calc_freematrix(S,@diff);
Kadv = calc_freematrix(S,@adv);


%%
[upcr,r] = random(upc);
kr = full(randomeval(k,r))
mat=FOUR_ISOT('k',kr,'b',[cx;cy],'stabilize','supg');
S = setmaterial(S,mat);
fr = bodyloadPG(S,[],'QN',FENODEFIELD(ones(nbnode,1)));
Kr = calc_rigi(S); 
ur = solve(Kr,fr);
figure(10)
clf
subplot(1,2,1)
plot(ur,S,'surface')
ax0=axis;
title('Petrov-Galerkin deterministe')
subplot(1,2,2)
plot(upcr,S,'surface')
axis(ax0)
title('Galerkin PC')


%%
uM = PCRADIALMATRIX([getnbddlfree(S),1],PC);

m=30;

stabil=0;

W=zeros(getnbddlfree(S),0);
for res=0:5
    Wadd=zeros(getnbddlfree(S),0);  
    l = rand(PC);
    l=l/norm(l);
    for i=1:min(m-getm(uM),8)
        kr = full(expect(k,l,l)/expect(l,l));  
        cr = [cx;cy];
        mat=FOUR_ISOT('k',kr,'b',cr,'stabilize','supg');
        S = setmaterial(S,mat);
        Kdiff = calc_freematrix(S,@diff);   
        Kadv  = calc_freematrix(S,@adv);  
        Kstab  = stabil*calc_freematrix(S,@stab);  
        Kr = Kdiff+Kadv+Kstab;
        if stabil
            fr = bodyloadPG(S,[],'QN',1)*expect(l)/expect(l,l);
        else
            fr = bodyload(S,[],'QN',1)*expect(l);    
        end
        fresr = fr;

        if getm(uM)>0
            fresr = fresr - Kdiff*expect(uM,k,l)/kr/expect(l,l);
            fresr = fresr - Kadv*expect(uM,l)/expect(l,l);
            fresr = fresr - Kstab*expect(uM,l)/expect(l,l);
        end
        w = solve(Kr,fresr);

        w=w/norm(w);
        w = w - Wadd*(Wadd'*w);
        fprintf('fun #%d, arnoldi residual = %d\n',i,norm(w));
        if norm(w)<1e-12
            break 
        end
        w=w/norm(w);
        Wadd=[Wadd,w];
        if getm(uM)>0
            respc = fpc-Kpc*uM;
        else
            respc = fpc; 
        end
        l = solve(w'*Kpc*w,w'*respc);
        l=l/norm(l);
    end

    W = [W,Wadd];
    Lam = solve(W'*Kpc*W,W'*fpc);
    uM =  PCRADIALMATRIX(W,size(uM),Lam);
    norm(Kpc*uM-fpc)/norm(fpc)
end




%%
uM = PCRADIALMATRIX([getnbddlfree(S),1],PC);
W=zeros(getnbddlfree(S),0);
m=4
stabil=1;

for i=1:m
    l0 = rand(PC); 
    l0 = l0/norm(l0);
    l=l0;
    for kk=1:6
        kr = full(expect(k,l,l)/expect(l,l));  
        cr = [cx;cy];
        mat=FOUR_ISOT('k',kr,'b',cr,'stabilize','supg');
        S = setmaterial(S,mat);
        Kdiff = calc_freematrix(S,@diff);   
        Kadv  = calc_freematrix(S,@adv);  
        Kstab  = stabil*calc_freematrix(S,@stab);  
        Kr = Kdiff+Kadv+Kstab;
        if stabil
            fr = bodyloadPG(S,[],'QN',1)*expect(l)/expect(l,l);
        else
            fr = bodyload(S,[],'QN',1)*expect(l)/expect(l,l);   
        end
        fresr = fr;
        if getm(uM)>0
            fresr = fresr - Kdiff*expect(uM,k,l)/kr/expect(l,l);
            fresr = fresr - Kadv*expect(uM,l)/expect(l,l);
            fresr = fresr - Kstab*expect(uM,l)/expect(l,l);
        end
        w = solve(Kr,fresr);
        w=w/norm(w);
        if getm(uM)>0
            respc = fpc-Kpc*uM;
        else
            respc = fpc; 
        end
        l = solve(w'*Kpc*w,w'*respc);  
%l=l/norm(l);
        errstagn = abs(norm(l)-norm(l0))/norm(l);
        fprintf('power iter %d , stagnation = %d\n',kk,errstagn);
        if errstagn<1e-2
            break 
        end
        l0=l;
    end
    W=[W,w];    
    Lam = solve(W'*Kpc*W,W'*fpc);
    uM =  PCRADIALMATRIX(W,size(uM),Lam);
    errres = norm(Kpc*uM-fpc)/norm(fpc) ;
    fprintf('nbfun %d , erreur res = %d\n',getm(uM),errres);

end


%%   
r = random(RANDVARS(PC));
%r = -1
uMr = randomeval(uM,r);
kr = full(randomeval(k,r));
upcr = full(randomeval(upc,r));
mat=FOUR_ISOT('k',kr,'b',[cx;cy],'stabilize','supg');
S = setmaterial(S,mat);
fr = bodyloadPG(S,[],'QN',1);
Kdiff = calc_freematrix(S,@diff);   
Kadv  = calc_freematrix(S,@adv);  
Kstab  = calc_freematrix(S,@stab);  
Kr = Kdiff+Kadv+Kstab;

ur = solve(Kr,fr);
figure(10)
clf
subplot(1,3,1)
plot(ur,S,'surface')
ax0=axis;
title('Petrov-Galerkin deterministe')
subplot(1,3,2)
plot(uMr,S,'surface')
axis(ax0)
title(['Galerkin GSD, order ' num2str(getm(uM))])
subplot(1,3,3)
plot(upcr,S,'surface')
axis(ax0)
title(['Galerkin PC'])
fprintf('\n')
fprintf('Erreur GSD-GAL = %d\n',norm(ur-uMr)/norm(ur));
fprintf('Erreur PC- GAL = %d\n',norm(ur-upcr)/norm(ur));
