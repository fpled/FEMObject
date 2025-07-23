r=10;
nt = 20;
P = POINT([ 0,0 ; 1,0 ; 0,2 ; 1,2 ; 1,1 ; 2,1 ; 2,2 ]);
S1 = mesh(DOMAIN(2,P(1),P(4)),r,r*2);
S2 = mesh(DOMAIN(2,P(5),P(7)),r,r);
S = union(S1,S2);
S = concatgroupelem(S);


k=46;

close all
display_=1;
rho=7800;
C=460;
%mat=FOUR_ISOT('c',rho*C,'k',k);

Tf = 1;
%Tf = 3e4;
t = linspace(0,Tf,nt+1);
dt=Tf/nt ;
T = TIMEMODEL(0,Tf,nt);


mat=FOUR_ISOT('c',1,'k',1);
S = setmaterial(S,mat);

L1 = LINE(P(6),P(7));
L2 = LINE(P(1),P(3));
L3 = LINE(P(3),P(7));
L4 = LINE(P(6),P(5));
L5 = LINE(P(2),P(5));
L6 = LINE(P(1),P(2));

S=final(S);
S=addcl(S,L1,'T',0);
S=addcl(S,L2,'T',0);
%S=addcl(S,L3,'T',0);
S=addcl(S,L4,'T',0);
S=addcl(S,L5,'T',0);


K=calc_rigi(S);
M=calc_mass(S);

f1 = surfload(S,L6,'QN',1);
f2 = surfload(S,L3,'QN',1);
f3 = bodyload(S,[],'QN',1);
f=f1+f2+f3;

q=K\f;

q=unfreevector(S,q);


if display_

    figure(1)
    clf

    plot(full(q),S,'surface')
    title('solution stationnaire')
%     axis image
    cax0=caxis;
end

%% resolution
N = DGTIMESOLVER(T,0);
f1t = -f1*rampe(N,0,Tf);
f2t = -f2*rampe(N,0,Tf);
f3t = f3*one(N);

ft = f1t+f2t+f3t;
qt = dsolve(N,ft,M,K);
qtmax = max(max(qt));
qtmin = min(min(qt));

N = setevolparam(N,'zlim',[qtmin,qtmax],'caxis',[qtmin,qtmax]);

if display_
    figure(15)
    evol(N,qt,S,'surface')
end
T = gettimemodel(T);

%% radial
testrad=1;
if testrad
    uref=qt;
    T = gettimemodel(N);
    urad = TIMERADIALMATRIX(size(qt),T);
    actu=1;

    Md=12;
    ftup = ft ; 
    Mt = getMmatrix(T);
    Dt=getDmatrix(T);
    for i=1:Md
        l0 = one(T);
        fprintf('--- fonction %d -----\n ',i)
        for j=1:7
            l=l0;
            A = M *integratemtimes(l,diff(l)) + K*integratemtimes(l,l);
%l=double(l0);
%A = M * (l*Mt*Dt*l') + K*(l*Mt*l');
            b = integratemtimes(ftup,l);
%b = double(ftup)*Mt*l';
            V = A\b;
            V = V/norm(V);
            l = dsolve(N,V'*ftup,V'*M*V,V'*K*V);
            fprintf('iteration puissance : erreur %d\n',norm(l-l0)/norm(l));
            l0=l;
%figure(10)
%plot(l);
%pause
%figure(11)
%plot(V,S);
%pause
        end
        urad = urad+TIMERADIALMATRIX(V,size(qt),l);    
        if actu==0
            ftup = ftup - (M*V)*diff(l) -(K*V)*l;
        else
            L = getL(urad);
            A = MULTIMATRIX(sparse(numel(M),numel(L)^2),size(M),[numel(L),numel(L)]);
            b = MULTIMATRIX(sparse(size(M,1),numel(L)),[size(M,1),1],[numel(L),1]);
            for ii=1:numel(L)
                for jj=1:numel(L)
                    A{ii,jj} = M *integratemtimes(L(ii),diff(L(jj))) + K*integratemtimes(L(ii),L(jj));
                    b{ii} = integratemtimes(ft,L(ii));
                end
            end
            A = assembleblock(A);
            b = assembleblock(b);
            V = A\b;
            ftup = ft;
            V=MULTIMATRIX(V,size(qt));

            for ii=1:numel(L)
                ftup = ftup - (M*V{ii})*diff(L(ii)) -(K*V{ii})*L(ii);
            end

            urad = TIMERADIALMATRIX(V,size(qt),L);    


        end
%ftup = ftup - TIMEMATRIX((M*V)*(Dt*double(l)')' +(K*V)*(Mt*double(l)')',T,size(V));
        err(i)=norm(urad-uref)/norm(uref);

        fprintf('erreur : %d\n',err(i));
    end


end
