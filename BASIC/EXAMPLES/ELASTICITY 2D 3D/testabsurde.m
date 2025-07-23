P=POINT([0,0;1,0;1,1;0,1]);
n=51;
MAT = ELAS_ISOT('E',1,'NU',0.3,'RHO',1);
D1=LINE(P(2),P(3));
D2=LINE(P(1),P(4));
S=mesh(DOMAIN(2,[0,0],[1,1]),n,n,MAT);
S=final(S,'renum');
S=addcl(S,D2,'U',0);
K=calc_rigi(S);
f=surfload(S,D1,'FY',1);
q=K\f;
clf;plot(S+0.1*q,'EdgeColor','r')
eval_sol(S,q,P(3),{'UX','UY'})
s=calc_sigma(S,q,'smooth');
plot(s,S,'compo','SMXX')
clear Ke;

clear rep

rep{1} = 1 ;
rep{2} = 2:n-1;
rep{3} = n;
rep{4} = [n+1 : n : n*(n-1)];
rep{6} = [2*n : n : n*(n-1)];
rep{7} = n*(n-1)+1;
rep{8} = n*(n-1)+[2:n-1];
rep{9} = n*n;
repbord=[];
for i=[1:4,6:9]
    repbord = union(repbord,rep{i});
end
rep{5} = setdiff(1:n^2,repbord);

% for e=1:S.nbelem
% repddl=S.groupelem{1}.numddl(e,:);
% Ke{e}=full(K(repddl,repddl));
% end

% for i=1:9
%     figure(1)
%     clf
%     ee=zeros(S.nbelem,1);ee(rep{i})=1;clf;plot(FEELEMFIELD(ee),S);
%     pause
% end

% for i=1:9
% eref = rep{i}(1);
% repddl=S.groupelem{1}.numddl(eref,:);
% Keref=full(K(repddl,repddl));
% repe=[];
% for e=1:S.nbelem
% if norm(Ke{e}-Keref)>1e-12
%     repe=[repe,e];
% end
% end
% ee=zeros(S.nbelem,1);ee(repe)=1;clf;plot(FEELEMFIELD(ee),S);
% pause
% end
% 
% q=MYDOUBLE(zeros(S.nbddl,1));
% m=length(S.ddlfree);
% rep=find(nonzeros(double(K)));
% A=sparse(size(K,1),size(K,2));
% A(rep)=double(K(rep)*eps);
% q(S.ddlfree)=((1+eps)*K(S.ddlfree,S.ddlfree))\f(S.ddlfree);
% figure(1);clf;plot(S+SPACEFUN(0.1*q),'EdgeColor','m')
