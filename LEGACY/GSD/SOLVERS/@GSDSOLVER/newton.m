function [u,result] = newton(GSD,L,b,calcAu,calcAtang,u0,Kt,solver)
% function u = newton(GSD,L,b,calcAu,calcAtang,u0,Kt,solver)
% solveur Newton GSD : resolution de A(u)=b
% GSD : GSDSOLVER
% L : NEWTONSOLVER
% b : second membre
% calcAu : fonction telle que calcAu(u) = A(u)
% calcAtang : fonction telle que calcAtang(u) = A'(u) (matrice tangente en u)
% u0 : initialisation
% Kt : matrice forc�e dans Newton modifiee
% solver : fonction donnant le solveur des systemes lineaires
%          solver(K,f) resout Ku=f

clock0 = clock;

type = getparam(L,'type');
fprintf('\n ---------------------------------------');
fprintf('\n ------------ Newton solver ------------')
fprintf('\n ---------------------------------------\n');

calcAu=fcnchk(calcAu);

try
    calcAtang=fcnchk(calcAtang);
catch
    warning('pas de fonction donnant la matrice tangente')
end

tol = getparam(L,'tol');
tolstagn = getparam(L,'tolstagn');
tolreact = getparam(L,'tolreact');
maxiter = getparam(L,'maxiter');
increment = getparam(L,'increment');


if nargin<6 || isempty(u0)
    u0 = zeros(size(b));
elseif isa(u0,'PCMATRIX')
    u0 = spectral_decomposition(u0,'tol',tol);
end

if nargin<7 || isempty(Kt) || ~strcmp(type,'manual')
    Kt = calcAtang(u0);
end


r0 = b-calcAu(u0);
err=norm(r0)/norm(b);
result.errorini=err;

if err<tol
    fprintf('Newton converged at initialization with residual error %.2d\n',err);
    return
end

fprintf('\nNEWTON initialization : residual = %d',err);
errstagn = 1;


for i=1:maxiter
    
    if increment
        if getparam(GSD,'reuse')
            [du,result.solver{i}] = solve(GSD,Kt,r0,u0);
        else
            [du,result.solver{i}] = solve(GSD,Kt,r0);
        end
        
        u = u0+du;
    else
        GSD = setparam(GSD,'tol',tol/10);
        GSD = setparam(GSD,'tolini',tol/10);
        if getparam(GSD,'reuse')
            [u,result.solver{i}] = solve(GSD,Kt,r0+Kt*u0,u0);
        else
            [u,result.solver{i}] = solve(GSD,Kt,r0+Kt*u0);
        end
        
    end
    
    %u = spectral_decomposition(u,'tol',tol/10,'nbfoncmax',length(u));
    
    %u=uniqueV(u);
    
    
    r = b-calcAu(u);
    
    
    err=norm(r)/norm(b);
    result.error(i)=err;
    if err<tol
        break
    else
        fprintf('\nNEWTON iteration #%3.d  : residual = %.2d',i,err);
    end
    
    errstagn = norm(u0-u)/norm(u+u0);
    
    if (strcmp(type,'full') || errstagn<tolstagn || err>tolreact) && ~strcmp(type,'manual')
        fprintf(' (actualisation tangent)')
        Kt = calcAtang(u);
    end
    
    u0 = u;
    r0 = r;
end

if err<tol
    fprintf('\nNewton converged at iteration #%d with residual error %.2d\n',i,err);
else
    fprintf('\nNewton stopped at iteration #%d with residual error %.2d\n',i,err);
end

added=[];for i=1:length(result.solver);added=[added,result.solver{i}.addedfunctions];end

fprintf('\n %3d added radial functions\n',sum(added))
disp(added)

fprintf('Elapsed time %f s\n',etime(clock,clock0))
