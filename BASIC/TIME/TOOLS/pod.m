function [u,result] = pod(b,varargin)
% function u = pod(b,varargin)
% Decomposition spectrale classique
% b : TIMEMATRIX ou TIMERADIALMATRIX
% u : TIMERADIALMATRIX
% u = sum_i U_i l_i
%     U_i vecteur, l_i fonctions temporelles
%     minimisant norme(b-sum_i U_i l_i)

paramradial.nbfoncmax = getcharin('nbfoncmax',varargin,50);
paramradial.tol = getcharin('tol',varargin,eps);
paramradial.pfixtol = getcharin('pfixtol',varargin,1e-3);
paramradial.pfixmax = getcharin('pfixmax',varargin,4);
uref = getcharin('reference',varargin);
display_ = ischarin('display',varargin);

if ~display_
    %    fprintf(' -> Decomposition spectrale ... ',class(b))
else
    fprintf('\n --- Decomposition spectrale d''un %s (POD)--- \n',class(b))
end

if istimeradial(b) && paramradial.nbfoncmax>getm(b)
    paramradial.nbfoncmax = getm(b);
end


errorpf=zeros(1,paramradial.pfixmax);
result.error=zeros(1,paramradial.nbfoncmax);
result.nfonctions = 0 ;
result.rayg=cell(1,paramradial.nbfoncmax);


T = gettimemodel(b);
sb = size(b);
l0 = zeros(1,T);
l = zeros(1,T);
U = zeros(sb(1),sb(2));
U0 = zeros(sb(1),sb(2));
radalpha = [];

if ischarin('metric',varargin)
    A = getcharin('metric',varargin);
else
    A = 1;
end


u = TIMERADIALMATRIX(sb,T);
bu=b ;

for j=1:paramradial.nbfoncmax
    
    l0 = rand(T);
    alpha0=norm(l0);l0=l0/alpha0;
    
    for i=1:paramradial.pfixmax
        
        U=integratemtimes(bu,l0);
        alpha = sqrt(U(:)'*A*U(:));
        U = U /alpha ;
        
        l = cell2mat(expand(U(:)'*A*bu(:)));
        alpha = norm(l);
        
        result.rayg{j}(i) = alpha;
        
        errorpf(i)=abs((result.rayg{j}(i)-(i>1)*result.rayg{j}(i-(i>1)))/result.rayg{j}(i));
        
        %errorpf(i)=double(sqrt(...
        %  (alpha0^2+alpha^2-2*alpha0*alpha*(U0(:)'*U(:))*(double(l0)*double(l)'))/...
        %                (alpha^2)));
        
        if display_
            fprintf('iteration %3d -> erreur = %.3e  \r',i,errorpf(i));
        end
        
        l0=l;U0=U;alpha0=alpha;
        if errorpf(i)<paramradial.pfixtol
            break
        end
    end
    
    
    
    radalpha=[radalpha,alpha];
    
    errorrayg = sqrt(alpha.^2/sum(radalpha.^2));
    
    if errorrayg< 10*eps ;
        fprintf('exacte avec %d fonctions -> erreur = %.3e \n',result.nfonctions,result.error(j))
        break
    end
    
    result.nfonctions = result.nfonctions + 1  ;
    
    u = u + TIMERADIALMATRIX(U,size(U),l) ;
    bu = bu - TIMERADIALMATRIX(U,size(U),l);
    
    
    if isempty(uref)
        result.error(j)=errorrayg   ;
    else
        result.error(j)=norm(expand(u-uref),A)/norm(expand(uref),A);
    end
    
    if display_
        fprintf('%d fonctions -> erreur = %.3e \n',result.nfonctions,result.error(j))
    end
    
    if result.error(j)<paramradial.tol
        break
    end
    
end


if ~display_
    %fprintf('nb fonctions %3d - erreur  : %3d \n',result.nfonctions,result.error(j))
end
