function [approx,V,indxi,varargout] = eim(snapshots,opts)
% function [approx,V,indxi,varargout] = eim(snapshots);
% function [approx,V,indxi,varargout] = eim(snapshots,opts);
%
% * Available options
%   opts.max_basis (= 20) Maximum dimension of the space used for the
%                         approximation
%
%   opts.tol (= 1e-12) Tolerance on the absolute error in
%                         $\max_j opts.metric(X(:,j))$
%
%   opts.orth (= true) Orthogonalize the basis of the space used for the
%                         approximation with respect to the inner product
%                         associated to opts.dot if defined, to the identity
%                         otherwise
%
%   opts.metric (= @(x) norm(x,Inf)) Error measure, it is the norm of interest
%
%   opts.dot (= []) Matrix used to compute the projection on the reduced
%                         subspace. If opts.dot is empty, a standard EIM is
%                         computed.

if nargin == 1
    opts = struct();
end
opts = default_opts(opts);
max_basis = opts.max_basis;
isproj = ~isempty(opts.dot);
if isproj
    A = opts.dot;
end

% m computes the norm of each column of snapshots according to opts.metric
m = @(y) cellfun(@(x) opts.metric(x), num2cell(y,1));

ind_xi = zeros(max_basis,1);
if ~isproj
    ind_x = zeros(max_basis,1);
end

[norm_s, indxi(1)] = max(m(snapshots));
V = snapshots(:,indxi(1));
if ~isproj
    V = V/opts.metric(V);
    [~,indx(1)] = max(abs(V));
else
    V = V/sqrt(V'*A*V);
end


for l = 2:max_basis
    if ~isproj
        Q = V(indx(1:l-1),:);
        c = Q\snapshots(indx(1:l-1),:);
        approx = V*c;
    else
        approx = V*((V'*A*V)\(V'*A*snapshots));
    end
    r = snapshots-approx;
    [norm_r, indxi(l)] = max(m(r));
    if norm_r < opts.tol
        indxi(l)=[];
        break
    end
    if ~isproj
        [~,indx(l)] = max(abs(r(:,indxi(l))));
    end
    v = r(:,indxi(l));
    v = v / norm(v);
    if opts.orth
        % orthogonalization of v w.r.t. V
        if ~isproj
            v = v - V*(V'*v);
            v = v / norm(v);
        else
            v = v - V*(V'*A*v);
            v = v / sqrt(v'*A*v);
        end
    end
    V = [V v];
end

% update the approximation after the last iteration
if ~isproj
    approx = V*(V(indx,:)\snapshots(indx,:));
else
    approx = V*((V'*A*V)\(V'*A*snapshots));
end

if ~isproj
    varargout{1} = indx;
end

end

function opts = default_opts(opts)
    if ~isfield(opts,'max_basis'); opts.max_basis = 20; end
    if ~isfield(opts,'tol'); opts.tol = 1e-12; end
    if ~isfield(opts,'orth'); opts.orth = true; end
    if ~isfield(opts,'metric'); opts.metric = @(x) norm(x,Inf); end
    if ~isfield(opts,'dot'); opts.dot = []; end
end
