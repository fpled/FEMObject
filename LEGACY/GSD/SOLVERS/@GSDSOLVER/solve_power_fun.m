function [u,result] = solve_power_fun(GSD,n,PC,fM,FM,f,F,res,ureuse,varargin)
% function u = solve(GSD,fM,FM,f,F,res,v,varargin)
% 
% v : pour la reutilisation (MULTIMATRIX ou PCRADIALMATRIX ou double)

display_ = getparam(GSD,'display');

fprintf('\n GSD solveur ... ')
if display_
    fprintf('\n')
end

paramradial.nbfoncmax = getparam(GSD,'nbfoncmax');
paramradial.tol = getparam(GSD,'tol');
paramradial.pfixtol = getparam(GSD,'pfixtol');
paramradial.pfixmax = getparam(GSD,'pfixmax');
paramradial.reuse = getparam(GSD,'reuse');
paramradial.update = getparam(GSD,'update');
paramradial.orthoduringpfix = getparam(GSD,'orthoduringpfix');

errorindicator = getparam(GSD,'errorindicator');
uref = getcharin('reference',varargin);

if ~isempty(uref)
    errorindicator='reference';
end

P = getP(PC);

errorpf=zeros(1,paramradial.pfixmax);
%result.error=zeros(1,paramradial.nbfoncmax);
%result.errorres=zeros(1,paramradial.nbfoncmax);
%result.errorrayg=zeros(1,paramradial.nbfoncmax);
%result.errorraygi=zeros(1,paramradial.nbfoncmax);
%result.time = zeros(1,paramradial.nbfoncmax);
%result.rayg=cell(1,paramradial.nbfoncmax);
%result.Rayg=cell(1,paramradial.nbfoncmax);
result.nfonctions = 0 ;
result.addedfunctions = 0;


l0 = zeros(1,PC);
l = zeros(1,PC);
U = sparse(n,1);
stoU = sparse(n,0);
U0 = sparse(n,1);
radalpha = [];


%% D�finition des normes 

b = res([]);    
mynorm=@(b) norm(b);
normb = mynorm(b);


%% D�finition des solveurs
toliter = getparam(GSD,'toliter');
directsolve = getparam(GSD,'direct');

clock0 = clock;


%% R�utilisation
if (paramradial.reuse && nargin>=4 && ~isempty(ureuse))
error('pas programme')
else
flagreuse = true;
end


if flagreuse
   u = PCRADIALMATRIX([n,1],PC);
   bini = b;
   result.raygini=0;
   result.errorini=1;    
else
result.raygini=resultini.rayg;
result.Raygini=resultini.Rayg;

switch errorindicator
    case 'reference'
result.errorini=mynorm(u-uref)/mynorm(uref);
    case 'residual'
result.errorini=mynorm(bini)/normb;            
    case 'none'
result.errorini=1;                    
end
if display_
fprintf('  Initialisation : %3d fonctions',length(u))    
end
if ~strcmp(errorindicator,'rayleigh')
if display_
fprintf(' - erreur  : %3d \n',result.errorini)
end
if result.errorini<getparam(GSD,'tolini')
fprintf('  Convergence apres initialisation - erreur : %3d\n',result.errorini)
return
end
end  
end
    
result.nfonctions = length(u);

result.totaltime=etime(clock,clock0);


%% Construction des couples de fonctions
nbfoncmax = paramradial.nbfoncmax;
bu=b;
for j=1:nbfoncmax %%%%%%%%%%%%%%%% BOUCLE NOUVELLES FONCTIONS

result.nfonctions = result.nfonctions + 1  ;

switch getparam(GSD,'inittype')
    case 'one'
l0 = one(PC);
    case 'random'         
l0 = rand(PC);
    case 'allone'
l0 = ones(length(PC),1);   
l0 = PCMATRIX(l0,[1,1],PC);

end

alpha0=norm(l0);l0=l0/alpha0;

for i=1:paramradial.pfixmax %%%%%%%%%%%%%%%% BOUCLE POINT FIXE

U = FM(l0);

%rayg(i)=U'*f0;
if paramradial.orthoduringpfix
V = double(getV(u));
if size(V,2)>=1
U = U - V*(V'*U);
end
end
alpha = sqrt(U'*U);
U = U /alpha ;
l = fM(U);
l = expand(l);

errorpf(i)=norm(l-l0)/(1/2*norm(l)+1/2*norm(l0));

if display_ 
    fprintf('  iteration %3d -> erreur = %.3e  \r',i,errorpf(i)); 
end

l0=l;U0=U;alpha0=alpha;

if errorpf(i)<paramradial.pfixtol 
    break 
end

end %%%%%%%%%%%%%%%% END BOUCLE POINT FIXE

u = u + U.*l;

if ~paramradial.update  %%%%%%%%%%%%% NO UPDATE OF STOCHASTIC FUNCTIONS
bu = res(u);

radalpha=[radalpha, sparse(size(radalpha,1),1);...
          sparse(1,size(radalpha,2)),norm(l)];

else   %%%%%%%%%%%%%% UPDATING STOCHASTIC FUNCTIONS
    
V = double(getV(u));
l = f(V);    
l = expand(l);

u = PCRADIALMATRIX(V,[n,1],l);
bu = res(u);

end %%%%%%%%%%%%% FIN UPDATING



% CALCUL ERREUR
switch errorindicator
    case 'residual'
result.error(j) = mynorm(bu)/normb;
    case 'reference'
result.error(j)=mynorm(u-uref)/mynorm(uref);
    case 'none'
result.error(j)=1;        
end
if ischarin('testerror',varargin)
%result.errorA(j)=mynorm(u-uref,A)/mynorm(uref,A);
result.errorres(j)=mynorm(bu)/normb;
end


if display_
fprintf('  %d fonctions -> erreur = %.3e \n',result.nfonctions,result.error(j))
if ischarin('testerror',varargin)
fprintf('  %d fonctions -> erreurres = %.3e \n',result.nfonctions,result.errorres(j))
end
end

result.addedfunctions = j;
result.time(j)=etime(clock,clock0);

if result.error(j)<paramradial.tol 
if strcmp(errorindicator,'rayleigh')
u = normsort(u);
u=truncate(u,length(u)-1);
result.nfonctions = length(u);
result.addedfunctions = j-1;  
end
hasconv=true;    
fprintf('  Convergence avec %3d fonctions -> erreur = %.3e \n',result.nfonctions,result.error(j))    
break
end

hasconv=false;

end %%%%%%%%%%%%%%%% END BOUCLE NOUVELLES FONCTIONS


result.totaltime=etime(clock,clock0);

%%%%% ITERATION SOUS ESPACE %%%%%%
if getparam(GSD,'subspaceiteration')>0

u0 = u;
l = getL(u0);
for kk=1:getparam(GSD,'subspaceiteration')
if display_
    fprintf('subspace iteration %d \n',kk)
end

V = F(l);
V = qr(V,0);
l = f(V);
u = PCRADIALMATRIX(V,[n,1],l);
    
errorpf(kk)=norm(u0-u)/norm(u0);
u0=u; 
end

bu = res(u);


% CALCUL ERREUR
switch errorindicator
    case 'residual'
result.error(j) = mynorm(bu)/normb;
    case 'reference'
result.error(j)=mynorm(u-uref)/mynorm(uref);
end
if ischarin('testerror',varargin)
%result.errorA(j)=mynorm(u-uref,A)/mynorm(uref,A);
result.errorres(j)=mynorm(bu)/normb;
end

if display_
fprintf('  %d fonctions -> erreur = %.3e \n',result.nfonctions,result.error(j))
if ischarin('testerror',varargin)
fprintf('  %d fonctions -> erreurres = %.3e \n',result.nfonctions,result.errorres(j))
end
end

end
%%%%% FIN ITERATION SOUS ESPACE %%%%%%



%% Finalisation : decomposition spectrale de la solution 

if getparam(GSD,'finalunique')
tolfactSD = getparam(GSD,'finaluniquefacttol');       
u=uniqueV(u,paramradial.tol*tolfactSD);    

elseif getparam(GSD,'finalSD')
 tolfactSD = getparam(GSD,'finalSDfacttol');   
 u = spectral_decomposition(u,'tol',paramradial.tol*tolfactSD,...
     'nbfoncmax',getm(u));  

end



%% Recalcul de la base stochastique finale  
if getparam(GSD,'finalupdate') %%%% FINAL UPDATING
U = double(getV(u));
l = f(U);

u = PCRADIALMATRIX(U,[n,1],l);
bu = res(u);
switch errorindicator
    case 'reference'
result.error(j)=mynorm(u-uref)/mynorm(uref);
    otherwise
result.error(j) = mynorm(bu)/normb;
end

disp('  Reactualisation')
fprintf('  nb fonctions %3d - erreur  : %3d \n',result.nfonctions,result.error(j))
end %%%% END FINAL UPDATING


%% Verdict sur la convergence
if ~hasconv
fprintf('  Non-convergence avec %3d fonctions -> erreur = %.3e \n',result.nfonctions,result.error(j))    
end

if display_
fprintf('Elapsed time is = %.3f seconds', etime(clock,clock0))
end
