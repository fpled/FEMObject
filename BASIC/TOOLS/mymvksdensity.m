function f = mymvksdensity(x,xi,varargin)
% function f = mymvksdensity(x,xi,'Bandwidth',bw)
% Computes a probability density estimate f of the sample (multivariate) data
% in the n-by-d matrix x, evaluated at the points in xi, using the required
% name-value pair argument value bw for the bandwidth.
% n is the number of points. d is the number of dimensions.
% xi is a d-column matrix that specifies the values where the density estimate is
% to be evaluated. bw is a d-element vector specifying the bandwidth of
% the  kernel smoothing window in each dimension or a scalar value for
% all dimensions. f is the vector of density value estimates.
% The estimation is based on a product Gaussian kernel function.

if nargin<=2 || isempty(varargin)
    error('Required argument value for Bandwidth')
else
    if ischarin('Bandwidth',varargin)
        bw = getcharin('Bandwidth',varargin);
    else
        error('Required argument value for Bandwidth')
    end
end

[n,d] = size(x);
% n = size(x,1); % number of data points
% d = size(x,2); % stochastic dimension
m = size(xi,1); % number of evaluation points

if isscalar(bw)
    bw = bw*ones(1,d);
elseif any(size(bw(:))~=[d,1])
    error(message('stats:mvksdensity:BadBandwidth',d));
end

cutoff = 4;

blocksize = 3e4;

if n*m<=blocksize
    % For small problems, compute kernel density estimate in one operation
    f = ones(n,m);
    for i = 1:d
        if verLessThan('matlab','9.1') % compatibility (<R2016b)
            z = bsxfun(@minus,xi(:,i)',x(:,i))/bw(i);
        else
            z = (xi(:,i)'-x(:,i))/bw(i);
        end
        % z = (repmat(xi(:,i)',n,1)-repmat(x(:,i),1,m))/bw(i);
        fi = normpdf(z);
        f = f.*fi;
    end
    f = mean(f);
else
    % For large problems, try more selective looping
    if isinf(cutoff)
        f = zeros(1,m);
        for i = 1:m
            ftemp = ones(n,1);
            for j = 1:d
                if verLessThan('matlab','9.1') % compatibility (<R2016b)
                    z = bsxfun(@minus,xi(i,j),x(:,j))./bw(j);
                else
                    z = (xi(i,j)-x(:,j))./bw(j);
                end
                % z = (repmat(xi(i,j),n,1)-x(:,j))./bw(j);
                fj = normpdf(z);
                ftemp = ftemp.*fj;
            end
            f(i) = mean(ftemp);
        end
    else
        weight = 1/n;
        halfwidth = cutoff*bw;
        index = (1:n)';
        f = zeros(1,m);
        for i = 1:m
            Idx = true(n,1);
            for j = 1:d
                dist = xi(i,j) - x(:,j);
                currentIdx = abs(dist) <= halfwidth(j);
                Idx = currentIdx & Idx; % pdf boundary
            end
            nearby = index(Idx);
            if ~isempty(nearby)
                ftemp = ones(length(nearby),1);
                for j = 1:d
                    if verLessThan('matlab','9.1') % compatibility (<R2016b)
                        z = bsxfun(@minus,xi(i,j),x(nearby,j))./bw(j);
                    else
                        z = (xi(i,j)-x(nearby,j))./bw(j);
                    end
                    % z = (repmat(xi(i,j),length(nearby),1)-x(nearby,j))./bw(j);
                    fj = normpdf(z);
                    ftemp = ftemp.*fj;
                end
                f(i) = weight*sum(ftemp);
            end
        end
    end
end
f = f(:)./prod(bw);

% For very small problems
% % if verLessThan('matlab','9.1') % compatibility (<R2016b)
% %     z = bsxfun(@rdivide,bsxfun(@minus,x',permute(repmat(xi',1,1,n),[1 3 2])),bw(:));
% % else
% %     z = (x'-permute(repmat(xi',1,1,n),[1 3 2]))./bw(:);
% % end
% z = (repmat(x',1,1,m)-permute(repmat(xi',1,1,n),[1 3 2]))./repmat(bw(:),1,n,m);
% f = reshape(prod(normpdf(z)),n,m);
% f = mean(f);
% f = f(:)./prod(bw);
