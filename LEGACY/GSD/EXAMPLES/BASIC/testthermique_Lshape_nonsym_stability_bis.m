% quelques conclusions : 
% - le calcul de w=F(lambda) tel que A(wl,w^*l)=B(w^*l) (Galerkin classique)
%    ok avec calcul discret stabilise
%   En notant a(u,v) = k a0(u,v) + a1(u,v) 
%   le probleme continu est a(u,v)=b(v)
%   et s(u,v) = int tau*Ladv(v)*L(u) = int tau*Ladv(v)*Ladv(u) 
%   on resout a(w,w^*)+s(w,w^*) = b(w^*)
%   c'est comme si on rajoutait une viscosite artificielle
%   pour le probleme considere, le second membre n'est pas modifie
%   par la stabilisation car int tau*Ladv(v)*f = 0 pour f constant
%
% avec la stabilisation simplifie pour grand peclet (stabtype 2), 
% la forme bilineaire de stabilisation est deterministe
%
% - le calcul de lambda = f(w) tel que A(wl,wl^*) = B(wl^*)
%   semble necessiter une stabilisation ?
%   ou simplement un meilleur calcul des integrales (autre definition 
%   des residus, lissage solution sur du degre 2, ...) ?
%   dans le system r(xi)*l = p
%   avec p=b(w) deterministe
%   et r(xi)=k(xi)a0(w,w) + a1(w,w)
%   ca donne un mauvais lambda
% en rajoutant a r(xi) la quantite s(w,w), ca marche (c'est juste un terme 
% deterministe mais assez grand). Bizarre, est-ce que c'est parce 
% que avec les elements finis, les termes a0(w,w) + a1(w,w) ne calculent 
% pas bien ce qu'on veut au niveau continu ????

% Il faudrait comparer les choses en 1D, en comparant une methode Galerkin
% spectrale avec Galerin-EF
% normalement la spectrale marche et comparer les r(xi) donnes par les deux
% approches
% les deux approches donnant a peu pres les memes w (a condition de stabiliser le 
% calcul deterministe EF)

%% modele stochastique
stabtype=1;
r=10;
p=4;

P = POINT([ 0,0 ; 1,0 ; 0,2 ; 1,2 ; 1,1 ; 2,1 ; 2,2 ]);
S1 = mesh(DOMAIN(2,P(1),P(4)),r,r*2);
S2 = mesh(DOMAIN(2,P(5),P(7)),r,r);
S = union(S1,S2);
S=concatgroupelem(S);
S=convertelem(S,'TRI3');
secdet=0;

%RV = RVUNIFORM(1/10,1/1);
RV = RVUNIFORM(1/1000,1/100)
X = PCMODEL(RV,'pcg','order',p);
PC = getPC(X);

k=X{1};
cx=0;
cy=1;


mat=FOUR_ISOT('k',k,'b',[cx;cy],'stabilize',stabtype);
S = setmaterial(S,mat);
S=final(S);
L1 = LINE(P(6),P(7));
L2 = LINE(P(1),P(3));
L3 = LINE(P(3),P(7));
L4 = LINE(P(6),P(5));
L5 = LINE(P(2),P(5));
L6 = LINE(P(1),P(2));



Kpcunfree = calc_rigipc(S,PC);   
% conditions aux limites
S=addcl(S,L1,'T',0);
S=addcl(S,L2,'T',0);
S=addcl(S,L3,'T',0);
S=addcl(S,L4,'T',0);
S=addcl(S,L5,'T',0);
S=addcl(S,L6,'T',0);

fnode = ones(S.nbnode,1); 
f=bodyloadPG(S,[],'QN',FENODEFIELD(fnode));
Kpc = calc_rigipc(S,PC);   
fpc = f*one(PC);
upc = solve(Kpc,fpc);

nbnode = S.nbnode;

mat=FOUR_ISOT('k',1,'b',[cx;cy]);
S = setmaterial(S,mat);
Kdiff1 = calc_freematrix(S,@diff);
Kadv = calc_freematrix(S,@adv);
Mas = calc_massgeom(S);

%% modele grossier
r=20;
S1 = mesh(DOMAIN(2,P(1),P(4)),r,r*2);
S2 = mesh(DOMAIN(2,P(5),P(7)),r,r);
Sg = union(S1,S2);
Sg=concatgroupelem(Sg);
Sg=convertelem(Sg,'TRI3');
mat=FOUR_ISOT('k',k,'b',[cx;cy]);
Sg = setmaterial(Sg,mat);
Sg=final(Sg);
Kgpcunfree = calc_rigipc(Sg,PC);   
Sg=addcl(Sg,L1,'T',0);
Sg=addcl(Sg,L2,'T',0);
Sg=addcl(Sg,L3,'T',0);
Sg=addcl(Sg,L4,'T',0);
Sg=addcl(Sg,L5,'T',0);
Sg=addcl(Sg,L6,'T',0);
fgnode = ones(Sg.nbnode,1); 
fg=bodyloadPG(Sg,[],'QN',FENODEFIELD(fgnode));
Kgpc = calc_rigipc(Sg,PC);   
fgpc = fg*one(PC);
ugpc = solve(Kgpc,fgpc);
Masg = calc_massgeom(Sg);

%%
stabtype=2;

[upcr,r] = random(upc);
kr = full(randomeval(k,r))
if exist('ugpc')
    ugpcr = full(randomeval(ugpc,r))
end
mat=FOUR_ISOT('k',kr,'b',[cx;cy],'stabilize',stabtype);
S = setmaterial(S,mat);
fr = bodyloadPG(S,[],'QN',FENODEFIELD(ones(nbnode,1)));
Kr = calc_rigi(S); 
ur = solve(Kr,fr);
figure(10)
clf
subplot(1,3,1)
plot(ur,S,'surface')
ax0=axis;
title('Petrov-Galerkin deterministe')
subplot(1,3,2)
plot(upcr,S,'surface')
%axis(ax0)
title('Galerkin PC')
subplot(1,3,3)
if exist('ugpc')
    plot(ugpcr,Sg,'surface')
%axis(ax0)
    title('Galerkin PC - fin')
end
fprintf('kr = %.1d,  Erreur PC- GAL = %d\n',kr,norm(ur-upcr)/norm(ur));

%%
GSD = GSDSOLVER('type','power','nbfoncmax',1,'display',true);
uM = solve(GSD,Kpc,fpc);


%%
uM = PCRADIALMATRIX([getnbddlfree(S),1],PC);

m=8;

stabil=1;

W=zeros(getnbddlfree(S),0);
for res=0:5
    Wadd=zeros(getnbddlfree(S),0);  
    l = rand(PC);
    l=l/norm(l);
    for i=1:min(m-getm(uM),8)
        kr = full(expect(k,l,l)/expect(l,l));  
        cr = [cx;cy];
        mat=FOUR_ISOT('k',kr,'b',cr,'stabilize',stabtype);
        S = setmaterial(S,mat);
        Kdiff = calc_freematrix(S,@diff);   
        Kadv  = calc_freematrix(S,@adv);  
        Kstab  = stabil*calc_freematrix(S,@stab);  
        Kr = Kdiff+Kadv+Kstab;

        fvolM = fnode;
        if getm(uM)>0
            fvolM = fvolM - Mas\(Kpcunfree*unfreevector(S,uM));
        end
        fvolMl = expect(fvolM,l)/expect(l,l);
        if stabil
            fresr = bodyloadPG(S,[],'QN',FENODEFIELD(fvolMl));
        else
            fresr = bodyload(S,[],'QN',FENODEFIELD(fvolMl));   
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
uMg = PCRADIALMATRIX([getnbddlfree(Sg),1],PC);
W=zeros(getnbddlfree(S),0);
Wg=zeros(getnbddlfree(Sg),0);
m=1
stabil=1;
stabtype=2;

update = 0;
fnodeM = fnode;
fgnodeM = fgnode;
for i=1:m
    l0 = rand(PC); 
    l0 = l0/norm(l0);
    l=l0;
    for kk=1:6
        kr = full(expect(k,l,l)/expect(l,l));  
        cr = [cx;cy];
        mat=FOUR_ISOT('k',kr,'b',cr,'stabilize',stabtype);
        S = setmaterial(S,mat);
        Kdiff = calc_freematrix(S,@diff);   
        Kadv  = calc_freematrix(S,@adv);  
        Kstab  = stabil*calc_freematrix(S,@stab);  
        Kr = Kdiff+Kadv+Kstab;

        fvolMl = expect(fnodeM,l)/expect(l,l);
        if stabil
            fresr = bodyloadPG(S,[],'QN',FENODEFIELD(fvolMl));
        else
            fresr = bodyload(S,[],'QN',FENODEFIELD(fvolMl));   
        end 
        w = solve(Kr,fresr);
        w=w/norm(w);
        wgunfree=transfer(S,Sg,w);
        wg = freevector(Sg,wgunfree);

        l = solve(wgunfree'*Kgpcunfree*wgunfree,wgunfree'*Masg*fgnodeM); 

%l=l/norm(l);
        errstagn = abs(norm(l)-norm(l0))/norm(l);
        fprintf('power iter %d , stagnation = %d\n',kk,errstagn);
        if errstagn<1e-2
            break 
        end
        l0=l;
    end

    if update
        W=[W,w];   
        Wg=[Wg,wg];
        Lam = solve(W'*Kpc*W,W'*fpc);
        Lamg = solve(Wg'*Kgpc*Wg,Wg'*fgpc);
        uM =  PCRADIALMATRIX(W,size(uM),Lamg);
        uMg = PCRADIALMATRIX(Wg,size(uMg),Lamg);

    else
        uM = uM + w*l;
        uMg = uMg + wg*l;

    end
    fgnodeM = cell2mat(fgnode - Masg\(Kgpcunfree*unfreevector(Sg,uMg)));
    fnodeM = transfer(Sg,S,fgnodeM);    

    errres = norm(Kgpc*uMg-fgpc)/norm(fgpc) ;
    fprintf('nbfun %d , erreur res = %d\n',getm(uM),errres);

end


%%
stovalr=[]
stovalMr=[]
stovalpcr=[]
stour=zeros(size(uM,1));
%%
stabtype=2;

r = random(RANDVARS(PC));
%r = -1
kr = full(randomeval(k,r));
mat=FOUR_ISOT('k',kr,'b',[cx;cy],'stabilize',stabtype);
S = setmaterial(S,mat);
fr = bodyloadPG(S,[],'QN',1);
Kdiff = calc_freematrix(S,@diff);   
Kadv  = calc_freematrix(S,@adv);  
Kstab  = calc_freematrix(S,@stab);  
Kr = Kdiff+Kadv+Kstab;

ur = solve(Kr,fr); 
uMr = randomeval(uM,r);
upcr = full(randomeval(upc,r));
stour=[stour,ur];
stovalr=[stovalr,ur(100)]
stovalMr=[stovalr,uMr(100)]
stovalpcr=[stovalr,upcr(100)]

figure(10)
clf
subplot(1,3,1)
plot(ur,S,'surface')
axis(ax0)
%ax0=axis;
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
