
multi_inclusions
S = createddlnode(S,DDL('u'));
S = addcl(S,[],'u');

%%
fx = LINFORM(0,1);
fx = setselgroup(fx,getnbgroupelem(S));
fx = fx{S}(:);
out = fx;

Ax = cell(1,getnbgroupelem(S));
for i=1:getnbgroupelem(S)
    Ax{i} = BILINFORM(1,1);
    Ax{i} = setselgroup(Ax{i},i);
    Ax{i} = Ax{i}{S}(:,:);
end

RV = RANDVARS(RVUNIFORM(0.9,1.1),getnbgroupelem(S)-2);
p=10;
[X,PC] = PCTPMODEL(RV,'order',p,'pcg','typebase',2);
%[X,PC] = PCMODEL(RV,'order',p,'pcg','nomasse');
%

A = (Ax{1}+Ax{end})*one(PC);
for i=1:getnbgroupelem(S)-2
    A = A + Ax{i+1}*X{i}*i^(-4);
end
b = fx*one(PC);

%%

Asep = SEPMATRIX(A);
bsep = SEPVECTOR(b);

%  calcul progressif

PGD = SEPSOLVER(getdim(Asep),'maxorder',30,'tol',1e-6,'display',true,...
    'maxiter',10,'updatedim',2:getdim(Asep),'update',1);
[usep,resultref] = solve(Asep,bsep,PGD);
useppc = PCTPMATRIX(usep,PC,2:getdim(usep));
oseppc=out'*useppc;
osep = mtimes(out',usep,1);
oref=osep;
uref = usep;
urefpc = useppc;
orefpc = out'*useppc;
sobref = sobol_indices(orefpc);
figure(12)
clf
semilogy(resultref.error)

%save('save_test_separation_multinclusions','uref','urefpc','orefpc','sobref','oref')

%%
P = POLYCHAOS(getM(PC),6);
I = getindices(P);I=I(:,1:end-1);
n=size(usep);
ucoeff=zeros(size(I,1),1);
for k=1:size(I,1)
    phi=cell(1,getM(PC)+1);
    phi{1}=speye(n(1));
for nu=1:size(I,2)
   phi{nu+1} = zeros(1,n(nu+1));
   phi{nu+1}(I(k,nu)+1) = 1;
end
phi=SEPMATRIX(phi);
ucoeff(k)= norm(squeeze(mtimes(phi,usep)));
end

%
[~,num] = sort(ucoeff,1,'descend')
Isorted = I(num,:)
    

%%
dim = [1,8];
p=500;
color='b';
Ip = Isorted(1:p,dim);
figure(1)
clf
plot(Ip(:,1),Ip(:,2),'LineStyle','none',...
    'Marker','s','MarkerFaceColor',color,...
    'MarkerEdgeColor',color,'MarkerSize',16)
set(gca,'FontSize',14)
xlim([-0.9,7])
ylim([-0.9,7])
grid on
set(gca,'GridLineStyle',':')

name = 'dim18';
myprint('./',['poly_indices_multiindices_veryaniso_' name],{'jpeg','epsc2'})

