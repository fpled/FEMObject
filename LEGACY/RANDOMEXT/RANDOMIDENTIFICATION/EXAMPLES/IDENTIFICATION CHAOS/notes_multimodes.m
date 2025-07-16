fprintf('\n\n\n')
X = RVLOGNORMAL(1,0.3,0,'stat');
Y = RVNORMAL(2,0.3);
Xs = random(X,1,10000);
%Xs = [random(X,1,1000)+random(Y,1,1000)];

mylike = @(X,Xs) (likelihood(X,Xs));
fprintf('log likelihood reference = %.1f \n',mylike(X,Xs))

%
param = lognfit(Xs);
Xfit = RVLOGNORMAL(param(1),param(2));
fprintf('log likelihood fit = %.1f \n',mylike(Xfit,Xs))

%
for p=[3,5,10]
    for nn=[50]
        Xpc = eicdfpcproject(Xs,p,POLYHERMITE(),'nbgauss',nn);
        fprintf('log likelihood pc order %d = %.1f \n',p,mylike(Xpc,Xs))

%
        figure(10)
        clf
        pdfsampleplot(Xs,'bar','w')
        xlim([min(Xs),max(Xs)])
%
        hold on
        pdfplot(Xfit,'g','npts',60)
        pdfplot(Xpc,'r','nbs',1e5,'npts',60)
        pause(1)
    end
end

%%
mu1 = -0.5;mu2=0.5;
si1=.35;si2=.35
ns1=100000;
ns2=100000;
Xs = [random('normal',mu1,si1,1,ns1),random('normal',mu2,si2,1,ns2)];

for p=[3,7,15,30,50]
    Xpc = eicdfpcproject(Xs,p,POLYHERMITE(),'nbgauss',50);
    fprintf('log likelihood pc order %d = %.1f \n',p,mylike(Xpc,Xs))
%
    figure(10)
    clf
    pdfsampleplot(Xs,'bar','w')
%
    xlim([-2.5,2.5])
    hold on
    pdfplot(Xpc,'r','nbs',3e5,'npts',100,'axis',xlim,'ksdensity')
    pause(1)
end

%% conditional pdf
sep=0.3;
for p=[1,3,5,10]
    Xpcm = eicdfpcproject_multimodes(Xs,p,POLYHERMITE(),sep,'nbgauss',40);
    fprintf('log likelihood pc order %d = %.1f \n',p,mylike(Xpcm,Xs))
    figure(10)
    clf
    pdfsampleplot(Xs,'npts',20,'bar','w','EdgeColor','c')
%
    xlim([-2.5,2.5])
    hold on
    pdfplot(Xpcm,'r','nbs',3e5,'npts',100,'axis',xlim,'ksdensity')
    pause(1)
end

%% polyfe
x = linspace(min(Xs),max(Xs),10);
for p=[1,3,5,10]
    Xpcfe = eicdfpcproject(Xs,p,POLYFE(x),'nbgauss',30);
    fprintf('log likelihood pc fe order %d = %.1f \n',p,mylike(Xpcfe,Xs))
    figure(10)
    clf
    pdfsampleplot(Xs,'npts',20,'bar','w','EdgeColor','c')
%
    xlim([-2.5,2.5])
    hold on
    pdfplot(Xpcfe,'r','nbs',3e5,'npts',100,'axis',xlim,'ksdensity')
    pause(1)
end


%% decomposition d'une variable a priori sur 1 germe de dim 2

PC = POLYCHAOS(2,2);
a = allones(PC);
Xs = random(a,1,500);
%a = PCMATRIX([1,0.5,0.3,0.4,0.2,0.1],[1,1],PC);

for p=[2,3,4,5,8]
    Xpc = eicdfpcproject(Xs,p,POLYHERMITE(),'nbgauss',40);
    fprintf('log likelihood pc fe order %d = %.1f \n',p,mylike(Xpc,Xs))
    figure(2)
    clf
    pdfplot(a,'nbs',1e5,'k')
    hold on
    pdfplot(Xpc,'nbs',1e5,'b')
    pause(1)
end


