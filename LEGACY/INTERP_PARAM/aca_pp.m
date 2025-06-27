function [F,I,J,varargout] = aca_pp(f,M,N,opts)
% [F,I,J,vargout] = aca_pp(f,M,N,opts)
% [F,I,J,Iref,Jref,Fref] = aca_pp(f,M,N,opts)
%
% * Available options:
%   opts.max_basis (= 20) Maximum dimension of the space used for the
%   approximation
%
%   opts.tol (= 1e-12) Tolerance with respect to the error
%   criterion selected
%
%   opts.pivot_search (= 'max') Algorithm for the selection of the
%   pivots. Available methods are 'max' and 'random_row'
%
%   opts.err_crit (= 'test_set') Error criterion. Available options
%   are 'stagnation' and 'test_set'. The test_set options compute
%   the infinite error on random entries of the matrix.
%
%   opts.err_card_test_set (= M + N) Used with the 'test_set'
%   options, it approximately defines the size of the test set
%   (up to the number of duplicates).
%
%   opts.err_confidence_level (= 0.9999) Confidence level for the
%   bound of the error on based on the test set.
%


if nargin == 3
    opts = struct();
end
opts = default_opts(opts,M,N);
opts.max_rank = min([opts.max_rank,M,N]);

if strcmp(opts.err_crit,'test_set')
    n_ref = opts.err_card_test_set;

    indref = randi(M*N, n_ref,1);
    n_ref = numel(indref);
    [Iref, Jref] = ind2sub([M N], indref);

    fprintf('#test_set = %d\n', n_ref);
    fprintf('Confidence level = %.4f%%\n', 100*opts.err_confidence_level);

    Fref = zeros(n_ref,1);
    for k = 1:n_ref
        Fref(k) = f(Iref(k),Jref(k));
    end

    t_alpha = tinv(1-opts.err_confidence_level,n_ref-1);

    varargout{1} = Iref;
    varargout{2} = Jref;
    varargout{3} = Fref;
end

F = zeros(M,N);

Frow = @(i) eval_F_row(f,i,N);
Fcol = @(j) eval_F_col(f,j,M);

I = zeros(opts.max_rank,1);
J = I;

U = zeros(M, opts.max_rank);
V = zeros(N, opts.max_rank);

err_crit = opts.tol + 1;

r = 1;
while (r <= opts.max_rank) && (err_crit > opts.tol)
    switch opts.pivot_search
      case 'max'
        [v, i] = pivot_max_search(Frow, F, V, r, I);
      case 'random_row'
        [v, i] = pivot_random_row_search(Frow, F, I, M);
      otherwise
        error('not implemented')
    end
    [~,j] = max(abs(v));

    v = v/v(j);
    u = Fcol(j) - F(:,j);

    I(r) = i;
    J(r) = j;
    U(:,r) = u;
    V(:,r) = v;

    F = U(:,1:r)*V(:,1:r)';

    switch opts.err_crit
      case 'stagnation'
        err_crit = norm(u)*norm(v);
        fprintf('iter %d -- err_crit = %d\n',r,err_crit);
      case 'test_set'
        err2_sample = (F(indref)-Fref).^2;
        err2_mean = mean(err2_sample);
        err2_std = std(err2_sample);

        err2_bound = err2_mean-t_alpha*err2_std/sqrt(n_ref);

        err_mean = sqrt(M*N*err2_mean);
        err_crit = sqrt(M*N*err2_bound);

        fprintf('iter %d -- err_mean = %d -- err_bound %d\n',...
                r, err_mean, err_crit);
      otherwise
        error('not implemented')
    end
    r = r + 1;
end

I = I(1:r-1);
J = J(1:r-1);

end

function opts = default_opts(opts,M,N)
if ~isfield(opts,'max_rank'); opts.max_rank = 20; end
if ~isfield(opts,'tol'); opts.tol = 1e-12; end
if ~isfield(opts,'pivot_search'); opts.pivot_search = 'max'; end
% Error criterion
if ~isfield(opts,'err_crit'); opts.err_crit = 'stagnation'; end
if strcmp(opts.err_crit,'test_set')
    if ~isfield(opts,'err_card_test_set')
        opts.err_card_test_set = M + N;
    end
    if ~isfield(opts,'err_confidence_level')
        opts.err_confidence_level = 99.99/100;
    end
end
end

function R = eval_F_row(f, i, N)
    R = zeros(N,1);
    for j = 1:N
        R(j) = f(i,j);
    end
end

function C = eval_F_col(f, j, M)
    C = zeros(M,1);
    for i = 1:M
        C(i) = f(i,j);
    end
end

function [v, i] = pivot_random_row_search(Frow, F, I, M)
    i = I(1);
    k = 0;
    nv = 0;
    while (nv == 0) && (k <= 100);
        l = 0;
        while any(i==I) && (l <= 100)
            i = randi(M);
            l = l + 1;
        end
        if (l <= 100) && any(i == I)
            error('No pivot available')
        end
        v = Frow(i) - F(i,:)';
        nv = norm(v);
        k = k + 1;
    end
end

function [v, i] = pivot_max_search(Frow,F,V,r,I)
    if r == 1
        i = 1;
    else
        pivots_available = 1:size(F,1);
        pivots_available(I(1:r-1)) = [];
        [~,piv] = max(V(pivots_available,r-1));
        i = pivots_available(piv);
    end
    v = Frow(i) - F(i,:)';
end
