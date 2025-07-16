%% definition des lois bimodales
nbmodes=2;
x1 = setnumber(RVUNIFORM(0,nbmodes),1);
switch nbmodes
case 2 
    mu = [1,2];
    si = [0.2,0.4];
    mufun = @(x) (x>=0 & x<1).*mu(1)  +  (x>=1 & x<=2).*mu(2);
    sifun = @(x) (x>=0 & x<1).*si(1)  +  (x>=1 & x<=2).*si(2);
case 3
    mu = [0,1,2];
    si = [0.3,0.4,0.2];
    mufun = @(x) (x>=0 & x<1).*mu(1)  +  (x>=1 & x<=2).*mu(2) + (x>=2 & x<=3).*mu(3);
    sifun = @(x) (x>=0 & x<1).*si(1)  +  (x>=1 & x<=2).*si(2) + (x>=2 & x<=3).*si(3);
case 4
    mu = [-1,1,2,5];
    si = [0.2,0.4,0.1,0.2];
    mufun = @(x) (x>=0 & x<1).*mu(1)  +  (x>=1 & x<=2).*mu(2) + (x>=2 & x<=3).*mu(3) + (x>=3 & x<=4).*mu(4);
    sifun = @(x) (x>=0 & x<1).*si(1)  +  (x>=1 & x<=2).*si(2) + (x>=2 & x<=3).*si(3) + (x>=3 & x<=4).*si(4);
end

x2 = setnumber(CONDRANDVAR(@RVNORMAL,x1,@(x) mufun(x),@(x) sifun(x)),2);        

figure(10)
clf
pdfplot(x2,'npts',100,'nbs',1e5);

Xs = random(x2,1000,1)';
hold on
pdfsampleplot(Xs,'b','npts',70,'bar')

%%
figure(11)
[C,x]=cdfsampleplot(Xs);
x = x(1:3:end);
C = C(1:3:end);
clf
plot(x,C)
title('CDF')
DC = diff(C)./diff(x);
xm = (x(1:end-1)+x(2:end))/2;
DDC = diff(DC)./diff(xm);
xmm = (xm(1:end-1)+xm(2:end))/2;
figure(12)
title('CDF')
plot(xm,DC,'r')
title('DERIVEE SECONDE DE CDF dont les zeros sont les points de stationarite de la pdf')

%%
sep = [1.5];

[P,x,Xblock,mublock] = separatesamples(Xs,sep,'display');


%%
figure(50)
clf
pdfsampleplot(Xs,'b','npts',110,'bar')
hold on
pc = POLYCHAOS(RANDPOLYS(POLYFE(x),POLYHERMITE()),[0,2],'typebase',2);
L = Inf;
scanL = [];
%%
for i=1:10
    [Xpc0,L0] = PCidentification(pc,Xs,'randomsearch','nbestim',2e3,'rs',20,'fmin',15,'parametrized','nocenter');
    scanL = [scanL;L0];
    if L0<L
        Xpc = Xpc0;
        L=L0;
    end
    figure(50)
    pdfplot(Xpc0,getcourbestyles(i,'nomarker'),'nbs',1e5,'npts',100);
    pause(2)

end
%%
%%
for i=1:10
    [Xpc0,L0] = PCidentification(pc,Xs,'mublock',[],...
        'nbestim',5e3,'rs',15,'fmin',25,'parametrized','nocenter');
    scanL = [scanL;L0];
    if L0<L
        Xpc = Xpc0;
        L=L0;
    end
    figure(50)
    pdfplot(Xpc0,getcourbestyles(i,'nomarker'),'nbs',1e5,'npts',100);
    pause(2)
    hold on
end


%%
%
[Xpc,L0] = PCidentification(pc,Xs,'plotfun',20,'mublock',mublock,...
    'nbestim',3e4,'rs',10,'fmin',10,'parametrized','nocenter');

%%
figure(50)
hold on 
pdfplot(Xpc,'npts',100,'nbs',1e5,'r','LineWidth',3)


%%
Xpcm = eicdfpcproject_multimodes(Xs,2,POLYHERMITE(),[1.5],'nbgauss',30)
figure(50)
clf
pdfsampleplot(Xs,'b','npts',110,'bar')
hold on
pdfplot(Xpcm,'r','npts',100,'nbs',1e5,'LineWidth',3)

-likelihood(Xpcm,Xs,1e5)

[Xpcmupdate,L0] = PCidentification(getPC(Xpcm),Xs,...
    'nbsestim',1e5,'fmin',20,'nocenter','initialguess',Xpcm,'parametrized');
pdfplot(Xpcmupdate,'g','npts',100,'nbs',1e5,'LineWidth',3)


%%

Xs = random(x2,1000000,1)';
figure(50)
clf
%pdfsampleplot(Xs,'b','npts',200,'ksdensity')
%hold on
pdfsampleplot(Xs,'b','npts',120)
hold on
%%
Xs = random(x2,1000,1)';
figure(50)
clf
pdfsampleplot(Xs,'w','npts',100,'ksdensity','bar','FaceColor','w')
hold on
%pdfsampleplot(Xs,'b','npts',120)
%hold on
Xpc = eicdfpcproject(Xs,15,POLYHERMITE(),'nbgauss',30)
pdfplot(Xpc,'g','nbs',1e5,'npts',120)

Xpcm = eicdfpcproject_multimodes(Xs,3,POLYHERMITE(),[1.6],'nbgauss',30)
pdfplot(Xpcm,'m','nbs',1e5,'npts',120)

Xpc = eicdfpcproject(Xs,4,POLYFE(x),'nbgauss',20)
pdfplot(Xpc,'k','nbs',1e5,'npts',120)

