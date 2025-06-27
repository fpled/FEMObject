function [W] = solve_rank1approxinv(A,varargin)
% function [W] = solve_rank1approxinv(A)
% function [W] = solve_rank1approxinv(A,P)
% function [W] = solve_rank1approxinv(A,opts)
% function [W] = solve_rank1approxinv(A,P,opts)


if nargin == 1
    P = [];
    opts = [];
elseif nargin == 2
    if isa(varargin{1},'htensor');
        P = varargin{1};
        opts = struct();
    elseif isa(varargin{1},'struct');
        P = [];
        opts = varargin{1};
    else
        error('Wrong argument')
    end
elseif nargin == 3
    P = varargin{1};
    opts = varargin{2};
elseif nargin>3
    error('Wrong number of arguments');
end

order = ndims(A);
opts = defaultopts(opts,order);

dim2ind = A.dim2ind;
n = sqrt(size(A));

error_iter = zeros(1,opts.maxiter);

I = sepmatrixtohtensor(sepopereye(n));
if opts.residual
    B = transpose_htensor(A);
else
    B = I;
end
C = apply_mat_to_mat(A,B,n);

if ~isempty(P)
    B = B-apply_mat_to_mat(P,C,n);
end

Br = rank(B);
Br = Br(dim2ind);
BU = matricize_leaves(B);

Cr = rank(C);
Cr = Cr(dim2ind);
CU = matricize_leaves(C);


W = I;
a = W.B{1};

WCW = cell(order,1);
BW = cell(order,1);



% cellfun works only for dense matrices :(
WU = matricize_leaves(W);
for mu = 1:order
    WCW{mu} = zeros(Cr(mu),1);
    WWmu = WU{mu}{1}'*WU{mu}{1}; % because dot(WC,W) = dot(C,W'W)
    for k = 1:Cr(mu)
        WCW{mu}(k) = sum(sum(CU{mu}{k}.*WWmu));
    end
    BW{mu} = zeros(Br(mu),1);
    for k = 1:Br(mu)
        BW{mu}(k) = sum(sum(BU{mu}{k}.*WU{mu}{1}));
    end
end


for i = 1:opts.maxiter
    a0 = a;
    for lambda = 1:order
        fprintf('%d ',lambda);
%         % Reduction of the core tensor
%         alpha = C.a;
%         m = size(alpha);
%         for mu = setdiff(1:order,lambda)
%             alpha = ttv(alpha,WCW{mu},mu);
%             m(mu) = 1;
%             alpha = reshape(alpha,m);
%         end
%         alpha = alpha(:);
%         % Construction of the operator for the dimension lambda
%         Q = alpha(1)*C.u{lambda}{1};
%         for j = 2:C.r(lambda)
%             Q = Q + alpha(j)*C.u{lambda}{j};
%         end
        Q = C;
        for mu = setdiff(1:order,lambda)
            Q.U{dim2ind(mu)} = WCW{mu}';
        end
        %Q = reshape(full(Q),n(lambda),n(lambda));
        % Previous line is OK, but it is memory consuming
        Q = reduce(Q,lambda);
        Q = reshape(Q,n(lambda),n(lambda));
        
%         % Reduction of the core tensor
%         beta = B.a;
%         m = size(beta);
%         for mu = setdiff(1:order,lambda)
%             beta = ttv(beta,BW{mu},mu);
%             m(mu) = 1;
%             beta = reshape(beta,m);
%         end
%         beta = beta(:);
%         % Construction of the operator for the dimension lambda
%         Z = beta(1)*B.u{lambda}{1};
%         for j = 2:B.r(lambda)
%             Z = Z + beta(j)*B.u{lambda}{j};
%         end
        Z = B;
        for mu = setdiff(1:order,lambda)
            Z.U{dim2ind(mu)} = BW{mu}';
        end
        % Z = reshape(full(Z),n(lambda),n(lambda));
        % Previous line is OK, but it is memory consuming
        Z = reduce(Z,lambda);
        Z = reshape(Z,n(lambda),n(lambda));
        if opts.symmetric
            if opts.sparse(lambda)
                v = gspai(Z,Q,'power',5,0.1);
                v = (v+v')/2;
            else
                if size(Q,1) < 1500
                    v = lyapsym(Q,(Z+Z'));
                else
                    v = full(Z)/full(Q);
                    v = (v+v')/2;
                end
            end
        else
            if opts.sparse(lambda)
                v = gspai(Z,Q,'power',5,0.1);
            else
                v = full(Z)/full(Q);
            end
        end
        a = norm(v,'fro');
        v = v/a;
        %W.a = tensor(a,ones(1,order));
        W.U{dim2ind(lambda)} = v(:);
        WWlambda = v'*v;
        for k = 1:Cr(lambda)
            WCW{lambda}(k) = sum(sum(CU{lambda}{k}.*WWlambda));
        end
        for k = 1:Br(lambda)
            BW{lambda}(k) = sum(sum(BU{lambda}{k}.*v));
        end
    end
    fprintf('\n');
    W.B{1} = a;
    error_iter(i) = abs(a-a0)/abs(a+a0);
    fprintf('Iter %d -- error %d\n',i,error_iter(i));
    if error_iter(i) < opts.stagcrit
        break;
    end
end

if opts.display==1
    if error_iter(i) < opts.stagcrit
        fprintf('     CV    - ');
    else
        fprintf('    als no CV - ');
    end
    toprint=['iter #' num2str(i,'%02.f') ' - error #' num2str(error_iter(i),'%.2e') ];
    disp(toprint);
end


end

function opts = defaultopts(opts,order)
    if ~isfield(opts,'maxiter');opts.maxiter = 10;end
    if ~isfield(opts,'stagcrit');opts.stagcrit = 5e-2;end
    if ~isfield(opts,'residual');opts.residual = 0;end
    if ~isfield(opts,'display');opts.display = 1;end
    if ~isfield(opts,'symmetric');opts.symmetric = 0;end
    if ~isfield(opts,'sparse');opts.sparse = zeros(1,order);end
end
