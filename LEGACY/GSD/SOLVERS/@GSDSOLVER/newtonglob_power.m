function [u,result] = newtonglob(GSD,N,S,b,X,varargin)
% function [u,result] = newtonglob(GSD,N,S,b,X,varargin)

fprintf('\nGSD solveur ... \n')
uref = getcharin('reference',varargin);
paramradial.nbfoncmax = getparam(GSD,'nbfoncmax');
paramradial.tol = getparam(GSD,'tol');
paramradial.pfixtol = getparam(GSD,'pfixtol');
paramradial.pfixmax = getparam(GSD,'pfixmax');
errorindicator = getparam(GSD,'errorindicator');
toliter = getparam(GSD,'toliter');
directsolve = getparam(GSD,'direct');
display_ = getparam(GSD,'display');
update = getparam(GSD,'update');


if ~isempty(uref)
    errorindicator='reference';
end


PC = getPC(X);
P = getP(PC);
n=size(b,1);

errorpf=zeros(1,paramradial.pfixmax);
result.error=zeros(1,paramradial.nbfoncmax);
result.time = zeros(1,paramradial.nbfoncmax);
result.nfonctions = 0 ;
result.addedfunctions = 0;

l0 = zeros(1,PC);
l = zeros(1,PC);
U = sparse(n,1);
U0 = sparse(n,1);
radalpha = [];

nbfoncmax = paramradial.nbfoncmax;


clock0 = clock;


u = PCRADIALMATRIX([n,1],PC);

N = setparam(N,'increment',true)
N = setparam(N,'stopini',false);


for j=1:nbfoncmax %%%%%%%%%%%%%%%% BOUCLE NOUVELLES FONCTIONS

result.nfonctions = result.nfonctions + 1  ;

l0 = ones(length(PC),1);l0=PCMATRIX(l0,[1,1],PC);
alpha0=norm(l0);l0=l0/alpha0;
N = setparam(N,'tol',paramradial.pfixtol/10);

for i=1:paramradial.pfixmax %%%%%%%%%%%%%%%% BOUCLE POINT FIXE

f0=expect(b,l0);

U = solve(N,f0,@(U) expect(calc_fintpc(S,(u+U*l0),X),l0),...
    @(U) expect(calc_rigitangpc(S,(u+U*l0),X),l0,l0),zeros(n,1));

alpha = sqrt(U'*U);
U = U /alpha ;

fU = U'*b;

l = solve(N,fU,@(l) U'*calc_fintpc(S,(u+U*l),X),...
    @(l) U'*calc_rigitangpc(S,(u+U*l),X)*U,zeros(1,PC));
%,[],@(a,f) cgs(a,f,getparam(N,'tol')/100));

result.rayg{j}(i) = expect(l,PCMATRIX(fU));
alpha = norm(l);
l=l/alpha;
%errorpf(i)=double(sqrt(...
%  (alpha0^2+alpha^2-2*alpha0*alpha*(U0'*Am*U)*(double(l0)*double(l)'))/...
%                (alpha^2)));
errorpf(i)=abs((result.rayg{j}(i)-(i>1)*result.rayg{j}(i-(i>1)))/result.rayg{j}(i));
            
if display_ 
    fprintf('  iteration %3d -> erreur = %.3e  \r',i,errorpf(i)); 
end

l0=l;U0=U;alpha0=alpha;

if errorpf(i)<paramradial.pfixtol 
    break 
end

end %%%%%%%%%%%%%%%% END BOUCLE POINT FIXE

u = u + U.*(alpha*l);

if update & getm(u)>1
%N = setparam(N,'tol',((j==1)/10+(j>1)*result.error(j-(j>1)))/10);
N = setparam(N,'tol',paramradial.tol/10);
if display_
fprintf('  reactualisation')
end
V = double(getV(u));
fU = V'*b;
l0 = zeros(getm(u),X);
l = solve(N,fU,...
    @(l) V'*calc_fintpc(S,PCRADIALMATRIX(V,[n,1],expand(l)),X),...
    @(l) V'*calc_rigitangpc(S,PCRADIALMATRIX(V,[n,1],expand(l)),X)*V,l0);   

u = PCRADIALMATRIX(V,[n,1],l);
end

    
switch errorindicator
    case 'residual'
res = b - calc_fintpc(S,(u),X);
result.error(j)=norm(res)/norm(b);
    case 'reference'
result.error(j)=norm(u-uref)/norm(uref);
end

N = setparam(N,'tol',result.error(j)/10);

if display_
fprintf('  %d fonctions -> erreur = %.3e \n',result.nfonctions,result.error(j))
end

result.addedfunctions = j;
result.time(j)=etime(clock,clock0);
if result.error(j)<paramradial.tol 
fprintf('  Convergence avec %3d fonctions -> erreur = %.3e \n',result.nfonctions,result.error(j))    
return
end

end %%%%%%%%%%%%%%%% END BOUCLE NOUVELLES FONCTIONS


fprintf('  Non-convergence avec %3d fonctions -> erreur = %.3e \n',result.nfonctions,result.error(j))    

result.totaltime=etime(clock,clock0);


