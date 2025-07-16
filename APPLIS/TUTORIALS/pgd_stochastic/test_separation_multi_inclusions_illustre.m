% ILLUSTRATIONS for the diffusion
% problem with multiple inclusions with random properties
% test_separation_multi_inclusions.m
mypath = './';
%%

figure(4)
clf
pdfplot(RV{1},'b','LineWidth',2)
hold on
pdfplot(RV{2},'g','LineWidth',2)
set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca,'FontSize',14)
legend('\alpha_1','\alpha_i, i= 2...8')
myprint(mypath,'ex_inclusions_rvs',{'epsc2'})
%%
figure(1)
clf

h1=plot(S,'EdgeColor','none','selgroup',[1,10],'FaceColor','r');
h2=plot(S,'EdgeColor','none','selgroup',2,'FaceColor','b');
h3=plot(S,'EdgeColor','none','selgroup',3:9,'FaceColor','g');
%h3=plot(S,'EdgeColor','none','selgroup',10,'FaceColor','r','');
plot(create_boundary( keepgroupelem(S,10)),'EdgeColor','w');
myprint(mypath,'ex_inclusions_geometry',{'jpeg'})
plot(S)
myprint(mypath,'ex_inclusions_mesh',{'jpeg'})
plot(create_boundary( keepgroupelem(S,10)),'EdgeColor','w');

%%

figure(3)
clf
plot(mean(urefpc),S);
myprint(mypath,'ex_inclusions_mean',{'jpeg'})
figure(4)
clf
plot(std(urefpc),S);
myprint(mypath,'ex_inclusions_std',{'jpeg'})

%%
i=2;
figure(3)
clf
plot(sobol_indices(urefpc,i),S);
myprint(mypath,['ex_inclusions_sob' num2str(i)],{'jpeg'})


%% resolution PGD Galerkin
update = 10;
PGD = SEPSOLVER(getdim(Asep),'maxorder',30,'tol',1e-5,'display',true,...
    'maxiter',10,'updatedim',2:getdim(Asep),'update',update,'storeiter',1,'itercritupdate',1e-8);
[usep,result] = solve(Asep,bsep,PGD);
useppc = PCTPMATRIX(usep,PC,2:getdim(usep));
oseppc=out'*useppc;
osep = mtimes(out',usep,1);

%% convergence des indices de sobol
%result=res1;
rrr=[];
for kk=1:length(result.u)
    kk
    output = PCTPMATRIX(result.u{kk},PC,2:getdim(result.u{kk}));
    output = out'*output;
%rrr = [rrr;sobol_indices(truncate(tempsol,1:kk))];
    rrr = [rrr;sobol_indices(output)];
end
figure(45)
plot(rrr,'LineWidth',2)
set(gca,'FontSize',14)
leg = {};
for k=1:8
    leg{k} = ['S_' num2str(k)];
end
ylim([0,0.4])
legend(leg{:},'Location','northeastoutside')
xlabel('rank')
%myprint(mypath,['ex_inclusions_progressive_galerkin_noupdate_pgd_sobconv'],{'epsc2'})
myprint(mypath,['ex_inclusions_progressive_galerkin_10update_pgd_sobconv'],{'epsc2'})


%%
i=15;
figure(12)
clf
pdfplot(truncate(oseppc,1:i),'nbs',1e5,'b','LineWidth',2,'ksdensity')
hold on
pdfplot(orefpc,'nbs',1e5,'k','LineWidth',2,'ksdensity')
legend(['Rank ' num2str(i)],'Reference')
myprint(mypath,['ex_inclusions_progressive_galerkin_3update_pgd_pdfconv_rank' num2str(i)],{'epsc2'})
%%
i=40;
figure(12)
clf
pdfplot(truncate(oseppc,1:i),'nbs',1e5,'b','LineWidth',2,'ksdensity')
hold on
pdfplot(truncate(osepgoalpc,1:i),'nbs',1e5,'r','LineWidth',2,'ksdensity')
hold on
pdfplot(orefpc,'nbs',1e5,'k','LineWidth',2,'ksdensity')
legend(['Rank ' num2str(i)],['Goal - Rank ' num2str(i)],'Reference')

%%
Qestimref = full(expect(orefpc));

for i=1:15
    Qestim = full(expect(truncate(oseppc,1:i)));
    err{1}(i)= abs(Qestim-Qestimref)/abs(Qestimref);
    Qestim = full(expect(truncate(osepgoalpc,1:i)));
    err{2}(i)= abs(Qestim-Qestimref)/abs(Qestimref);
end

figure(23)
clf
semilogy(err{1})
hold on
semilogy(err{2},'r')


%% convergence in L2 norm

PGD = SEPSOLVER(getdim(Asep),'maxorder',20,'tol',1e-5,'display',true,...
    'maxiter',15,'itercrit',1e-4,'updatedim',2:getdim(Asep),'update',2,...
    'reference',uref);
[usep,result] = solve(Asep,bsep,PGD);
up = getparam(PGD,'update');
resup{up+1} = result;

%%
figure(12)
clf
leg = {};
for k=0:length(resup)-1
    semilogy(resup{k+1}.error,getcourbestyles(k+1),'LineWidth',2)
    hold on
    set(gca,'FontSize',14)
    xlabel('rank')
    ylabel('error')
    leg{end+1} = [num2str(k) ' updates'];
end
legend(leg{:})

myprint(mypath,'ex_inclusions_convL2_updates',{'epsc2','jpeg'})
