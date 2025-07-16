%% 8D POISSON PROBLEM
dim=8;
nelem=60;

D=DOMAIN(1);
S=mesh(D,nelem);
S=createddlnode(S,DDL('u'));
S=addcl(S,POINT([0;1]),'u');

K=BILINFORM(1,1,1);
K=K{S}(:,:);
M=BILINFORM(0,0,1);
M=M{S}(:,:);

A=SEPMATRIX(dim);
for d=1:dim
    ad=cell(1,d);
    for dd=1:dim
        if dd==d
            ad{dd}=K;
        else
            ad{dd}=M;
        end
    end
    A=A+SEPMATRIX(ad);
end

x=getcoord(getnode(S));
x=double((x<2/3) & (x>1/3));
xx=LINFORM(0,x,0);
xx=xx{S}(:);
b=cell(1,dim);
for d=1:dim
    b{d}=xx;
end
b=SEPMATRIX(b);

%% GALERKIN PGD
PGD = SEPSOLVER(A.dim,'maxorder',8,...
    'errorindicator','residual');
[u0,r0] = solve(A,b,PGD);

%% OPERATORS CONVERSION TO HTENSOR FORMAT
Ah = sepmatrixtohtensor(A);
Ah = sparse_leafs(orthog(full_leafs(Ah)));

%%
opts.max_rank = max(rank(Ah));
opts.rel_eps = 1e-10;
opts.disp_errtree = 1;
Ah = truncate(Ah,opts);

bh = sepmatrixtohtensor(b);
bh = orthog(bh);

%% GALERKIN HPGD
tic
[u1,r1] = solve_hpgd_prog(Ah,bh,PGD);
toc

%% GALERKIN HPGD WITH LEAVES UPDATE
HPGD = setparam(PGD,'updateleaves',1);
tic
[u2,r2] = solve_hpgd_prog(Ah,bh,HPGD);
toc
%%
figure;
semilogy(r0.error)
hold on
semilogy(r1.error,'r')
semilogy(r2.error,'g')

legend('PGD','HPGD','HPGD with leaves update',...
    'Location','southwest');
xlabel('Rank of the approximation')
ylabel('Relative Residual')

%% DIFFERENT TREE TEST
% see define_tree for predefined trees
% First dim separated
Ah = sepmatrixtohtensor(A,'first_separate');
Ah = sparse_leafs(orthog(full_leafs(Ah)));
opts.max_rank = max(rank(Ah));
opts.rel_eps = 1e-10;
opts.disp_errtree = 1;
Ah = truncate(Ah,opts);
bh = sepmatrixtohtensor(b,'first_separate');
bh = orthog(bh);
[u3,r3] = solve_hpgd_prog(Ah,bh,PGD);

% TT tree
Ah = sepmatrixtohtensor(A,'TT');
Ah = sparse_leafs(orthog(full_leafs(Ah)));
opts.max_rank = max(rank(Ah));
opts.rel_eps = 1e-10;
opts.disp_errtree = 1;
Ah = truncate(Ah,opts);

bh = sepmatrixtohtensor(b,'TT');
bh = orthog(bh);
[u4,r4] = solve_hpgd_prog(Ah,bh,PGD);

%%
figure
semilogy(r1.error)
hold on
semilogy(r3.error,'r')
semilogy(r4.error,'g')

legend('Balanced tree','First dim separated','TT',...
    'Location','southwest');
xlabel('Rank of the approximation')
ylabel('Relative Residual')

