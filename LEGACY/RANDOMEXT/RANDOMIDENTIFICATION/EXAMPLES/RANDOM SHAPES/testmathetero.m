L = 30;
rmin = 1;
rmax = 2;
rmean = rmax/2+rmin/2;
alpha=.5;
N = floor(alpha*L^2/pi/rmean^2);

tir = 2;
switch tir
case 1
    c = lhsdesign(N,2)*L;
    cx=c(:,1);
    cy=c(:,2);
    r = random(RVUNIFORM(rmin,rmax),N,1);
case 2
    RV = RANDVARS(RVUNIFORM(0,L),RVUNIFORM(0,L),RVUNIFORM(rmin,rmax))
    cx=random(RV{1},N,1);
    cy=random(RV{2},N,1);
    r =random(RV{3},N,1);
end

%
figure(1)
clf
plot(DOMAIN(2,[0,0],[L,L]));
for i=1:N
    plot(CIRCLE(cx(i),cy(i),r(i)))
end
axis square

%%
l = L;
nbox = 1;
cxbox = random(RVUNIFORM(l/2,L-l/2),nbox,1);
cybox = random(RVUNIFORM(l/2,L-l/2),nbox,1);
dlbox = [cxbox-l/2,cybox-l/2];
urbox = [cxbox+l/2,cybox+l/2];

figure(1)
for i=1:nbox
    plot(DOMAIN(2,dlbox(i,:),urbox(i,:)),'EdgeColor','r')
end
%%

[S,X,Y] = mesh(DOMAIN(2,[0,0],[l,l]),1000,1000);
ls = cell(1,nbox);
scanbox = 1:nbox;
lss = zeros(numel(X),nbox);
for i=1:nbox
    pourcentage(i,nbox,50);
    rep = find(cx<=urbox(i,1)+r & cy<=urbox(i,2)+r ...
        &cx>=dlbox(i,1)-r & cy>=dlbox(i,2)-r);
    j = 1;
    lsi = -lseval(LSCIRCLE(cx(rep(j))-dlbox(i,1),cy(rep(j))-dlbox(i,2),r(rep(j))),S);
    for j=2:length(rep)
        lstemp = -lseval(LSCIRCLE(cx(rep(j))-dlbox(i,1),cy(rep(j))-dlbox(i,2),r(rep(j))),S);
        lsi = intersect(lsi,lstemp); 
    end
    lss(:,i)=lsi;
end

%%

[X,Y] = meshgrid(linspace(0,L,1000),linspace(0,L,1000));
for k=1:length(cx)

    pourcentage(k,length(cx),1)
    lst =  r(k) - sqrt((X-cx(k)).^2+(Y-cy(k)).^2);
    if k==1
        lss=lst;
    else
        lss = max(lst,lss);  
    end

end

figure(44)
clf
surf(X,Y,double(lss<0),'EdgeColor','none')
clear lst

%%
A = double(lss<0);

figure(5)
clf
subplot(1,2,1)
surf(X,Y,A,'EdgeColor','none')
view(3)


s=svd(A);
%s=diag(D);
err=sqrt(1-cumsum(s.^2)/sum(s.^2));
figure(13)
semilogy(err,'b--')
hold on
clear A
clear V
clear D
clear U


coeff = 20;
fun = @(x) (x>0).*(1-1/2*exp(-coeff*abs(x)))+(x<=0).*(1/2*exp(-coeff*abs(x)));
Af=fun(lss);
figure(5)
subplot(1,2,2)
surf(X,Y,Af,'EdgeColor','none')
view(3)

s=svd(Af);
%s=diag(D);
errf=sqrt(1-cumsum(s.^2)/sum(s.^2));

rep = find(errf<1e-2);
rep(1)

figure(13)
semilogy(errf,'g--')
hold on


%%
figure(5)
for i=1:nbox
    clf
    plot(S)
    contourplot(LEVELSET(lss(:,i)),S,'Color','r')
    pause(.5)
end    

%%
A = reshape(lss(:,1),size(X));

[U,D,V]=svd(A);
s=diag(D);
err=sqrt(1-cumsum(s.^2)/sum(s.^2));

for m=1:10:size(D,1)
    figure(5)
    clf
    lssapprox = U(:,1:m)*D(1:m,1:m)*V(:,1:m)'; 
    plot(S)
    contourplot(LEVELSET(lssapprox(:)),S,'Color','r')
    title(['m = ' num2str(m) ' , error = ' num2str(err(m))])
    pause
end
%%

A = reshape(double(lss(:,1)),size(X));
figure(5)
clf
surf(X,Y,A,'EdgeColor','none')
coeff = 10;
fun = @(x) (x>0).*(1-1/2*exp(-coeff*abs(x)))+(x<=0).*(1/2*exp(-coeff*abs(x)));
%fun = @(x) x;
figure(5)
clf
surf(X,Y,fun(A),'EdgeColor','none')

A=fun(A);


[U,D,V]=svd(A);
s=diag(D);
err=sqrt(1-cumsum(s.^2)/sum(s.^2));

rep = find(err<1e-2);
rep(1)

figure(13)
semilogy(err,'r')
hold on
%%
for m=10:30:size(D,1)
    figure(5)
    clf
    lssapprox = U(:,1:m)*D(1:m,1:m)*V(:,1:m)'; 
    plot(S)
%contourplot(LEVELSET(lssapprox(:)),S,'Color','r')
    surf(X,Y,lssapprox,'EdgeColor','none')
    view(3)
    title(['m = ' num2str(m) ' , error = ' num2str(err(m))])
    pause(0.8)
end

%%
A = reshape(double(lss(:,1)<0),size(X));

[U,D,V]=svd(A);
s=diag(D);
err=sqrt(1-cumsum(s.^2)/sum(s.^2));


rep = find(err<1e-2);
rep(1)

figure(13)
semilogy(err,'b')
hold on
%%

for m=1:30:size(D,1)
    figure(5)
    clf
    lssapprox = U(:,1:m)*D(1:m,1:m)*V(:,1:m)'; 
    plot(S)
%contourplot(LEVELSET(lssapprox(:)),S,'Color','r')
    surf(X,Y,lssapprox,'EdgeColor','none')
    view(2)
    title(['m = ' num2str(m) ' , error = ' num2str(err(m))])
    pause(.5)
end
%%
lssm = mean(lss,2);
lsscentered = lss - repmat(lssm,1,size(lss,2));
[U,D,V] = svd(lsscentered,'econ');
%
figure(10)
clf
d = diag(D)
semilogy(d(1:end-1))
title('spectre operateur covariance')
