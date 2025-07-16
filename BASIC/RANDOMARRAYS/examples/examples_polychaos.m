
%% CHAOS POLYNOMIAL

% chaos d'hermite d'ordre 4 en dimension 1 
h1 = POLYHERMITE()
PC1 = POLYCHAOS(h1,4) % chaos d'ordre 2 et en dimension 1
getindices(PC1)

% chaos de legendre d'ordre 3 en dimension 2
h2 = RANDPOLYS(POLYLEGENDRE(),POLYLEGENDRE())  % polynomial multi-dimensionel
PC2 = POLYCHAOS(h2,3) % chaos d'ordre 2 et en dimension 1
getindices(PC2)


%% VARIABLES ALEATOIRES

help RVNORMAL 
a = RVNORMAL(1,2)
figure(1)
subplot(1,2,1)
pdfplot(a)
subplot(1,2,2)
cdfplot(a)

help RVUNIFORM
b = RVUNIFORM(-1,2)
figure(2)
subplot(1,2,1)
pdfplot(b)
subplot(1,2,2)
cdfplot(b)

help RVLOGNORMAL
c = RVLOGNORMAL(2,0.3,1,'stat')
figure(3)
subplot(1,2,1)
pdfplot(c)
subplot(1,2,2)
cdfplot(c)

help RVEXP
d = RVEXP(1,1)  % loi exponentielle, cas particulier de la loi Gamma 
figure(4)
subplot(1,2,1)
pdfplot(d)
subplot(1,2,2)
cdfplot(d)

%% decomposition d'une variable aleatoire sur le chaos 

apc = decomppc(a,'order',4);
full(apc)
figure(5)
clf
pdfplot(a)
hold on
pdfplot(apc,'r','nbs',1e5) % la pdf est calculer en faisant des tirages du developpement sur le chaos : on peut donner le nombre de tirages 'nbs' 

bpchermite = decomppc(b,'order',6)
bpclegendre = decomppc(b,POLYLEGENDRE(),'order',6)
figure(6)
clf
pdfplot(b)
hold on
pdfplot(bpchermite,'r','nbs',1e5)
pdfplot(bpclegendre,'g','nbs',1e5)


%% decomposition d'une fonction de variable aleatoire sur le chaos
% fonction f(a) = a^2-2 avec a decomposee sur le chaos
apc = decomppc(a,POLYHERMITE(),'order',4);
fun = @(a) a.^2-2;
fpc = decompfun(apc,[],[],fun) % le deuxieme argument est le nombre de points de gauss par dimension

%% creation d'un modele CHAOS POLYNOMIAL A PARTIR DE PLUSIEURS VARIABLES ALEATOIRES
X = PCMODEL(RANDVARS(a,d),'pcg','order',5) % chaos generalise
% on obtient les deux variables aleatoires a et d sur le chaos, notees X1
% et X2
figure(10)
clf
subplot(1,2,1)
pdfplot(X{1})
subplot(1,2,2)
pdfplot(X{2})
% decomposition de la fonction f(X1,X2) = [X1+2*X2;2*X1]  (un vecteur aleatoire)
fun = @(x) MULTIMATRIX([x(:,1)+2*x(:,2),2*x(:,1)]',[2,1]);
% utilisation de decompmatrix
help decompmatrix
fpc = decompmatrix(X,[],[],fun)
figure(11)
clf
subplot(1,2,1)
pdfplot(fpc(1))
subplot(1,2,2)
pdfplot(fpc(2))
figure(12)
multipdfplot(fpc,'nbs',1e4,'npts',50)


%% FONCTION DE VARIABLE ALEATOIRE
apc = decomppc(a,POLYHERMITE(),'order',4);
pc = getPC(apc);
r = RANDVARFUNCTION(@(x) x.^2-2 ,apc);
rpc = project(r,pc)


%% CALCUL DES POINTS DE GAUSS POUR UN CHAOS
PC = POLYCHAOS(POLYHERMITE(),3);
H = RANDPOLYS(PC)
gauss = calc_gausspoints(PC,8);
figure(100)
plot(gauss.coord,zeros(1,length(gauss.w)),'o')
gauss.w
grid on

PC = POLYCHAOS(RANDPOLYS(POLYHERMITE(),POLYLEGENDRE()),3);
H = RANDPOLYS(PC)
gauss = calc_gausspoints(PC,5);
figure(100)
plot(gauss.coord(:,1),gauss.coord(:,2),'o')
gauss.w
grid on


%% VARIABLES CONDITIONNELLES

H = setnumber(RVGUMBEL(28.5,1.9),1);
T = CONDRANDVAR(@RVNORMAL,H,@(H) sqrt(5.88*H+38.82),@(H) 0.2*sqrt(5*H));
RV = RANDVARS(H,T);
X = PCMODEL(RV,'order',4)

figure(1)
plot(H,'pdf');
hold on
plot(X{1},'pdf','r','nbs',1e4)

figure(2)
plot(T,'pdf');
hold on
plot(X{2},'pdf','r','nbs',1e4)

figure(3)
subplot(1,2,1)
plot(T,'multipdf','npts',20,'EdgeColor','none')
view(2)
ax = axis;
subplot(1,2,2)
plot(X(:),'multipdf','npts',20,'EdgeColor','none')
view(2)
axis(ax)
