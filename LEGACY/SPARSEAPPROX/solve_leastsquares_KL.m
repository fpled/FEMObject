function [u,err,r] = solve_leastsquares_KL(A,y,regul,param,cv,k,tolKL,cvKL,kKL,varargin)
% function [u,err,r] = solve_leastsquares_KL(A,y,regul,param,cv,k,tolKL,cvKL,kKL,varargin)
% Perform a piori empirical Karhunen-Loeve (KL) decomposition, solve (regularized) least-squares minimization problem
% and compute relative cross-validation error
% A : matrix containing the evaluations of basis functions
% y : response vector or matrix
% regul : type of regularization ('' or 'l0' or 'l1'), '' by default
% param : parameters for regularization, struct() by default
% cv : type of cross-validation procedure ('leaveout' or 'kfold'), 'leaveout' by default
% k : number of folds (only for k-fold cross-validation procedure), 10 by default
% tolKL : prescribed tolerance for truncated empirical KL decomposition
% cvKL : type cross-validation procedure for empirical KL decomposition ('leaveout' or 'kfold'), 'leaveout' by default
% kKL : number of folds (only for k-fold cross-validation), 10 by default
% 'ploterrorKL' in varargin : plot relative SVD truncation error w.r.t. SVD rank (true or false), false by default
%
% See also SPARSEAPPROX/solve_leastsquares

if nargin<3 || isempty(regul)
    regul = '';
end
if nargin<4 || isempty(param)
    param = struct();
end
if nargin<5 || isempty(cv)
    cv = 'leaveout';
end
if strcmp(cv,'kfold') && (nargin<6 || isempty(k))
    k = 10;
end
if nargin<7 || isempty(tolKL)
    tolKL = [];
end
if nargin<8 || isempty(cvKL)
    cvKL = 'leaveout';
end
if strcmp(cvKL,'kfold') && (nargin<9 || isempty(kKL))
    k = 10;
end

% Reduction step
N = size(y,1);
mu = mean(y,1);
% y_centered = y - repmat(mu,N,1);
y_centered = y - mu.*ones(N,1);
[U,S,V,errKL] = svdtruncate(full(y_centered'),tolKL);
r = length(diag(S));
% errKL = zeros(N,1);
% for k=1:N
%     y_centered_approx = U*S*V';
%     errKL(k) = norm(y_centered - y_centered_approx,'fro')/norm(y_centered,'fro');
% end
if ischarin('display',varargin)
    fprintf('\nSVD rank : r = %d\r',r);
end

% Plot relative SVD truncation error
if ischarin('ploterrorKL',varargin)
    plot_SVD_error(errKL);
end

% Least-squares minimization step
v = solve_leastsquares(A,V,regul,param,cv,k,varargin{:});

% Reconstruction step
u = U*S*v';
u(:,1) = u(:,1) + mu';

% Cross-validation step
if nargout>1
    fun = @(y,A) solve_leastsquares(A,y,regul,param,cv,k,varargin{:});
    err = calc_cv_error(y,A,fun,cvKL,kKL,varargin{:});
end

end

function plot_SVD_error(err,varargin)

p = ImprovedInputParser;
addParameter(p,'legend',true,@islogical);
addParameter(p,'grid',true,@islogical);
addParameter(p,'box',true,@islogical);
addParameter(p,'FontSize',16,@isscalar);
addParameter(p,'LineWidth',1,@isscalar);
parse(p,varargin{:})

figure('Name','Evolution of relative SVD truncation error indicator w.r.t SVD rank')
% set(gcf,'Name','Evolution of relative SVD truncation error indicator w.r.t SVD rank')
semilogy(1:numel(err),err,'-+k','LineWidth',p.Results.LineWidth);
if p.Results.grid
    grid on
end
if p.Results.box
    box on
end
set(gca,'FontSize',p.Results.FontSize)
xlabel('SVD rank')
ylabel('SVD truncation error')
if p.Results.legend
    legend('SVD rank versus SVD truncation error')
end

end

