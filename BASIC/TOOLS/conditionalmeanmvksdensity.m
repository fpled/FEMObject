function [my,bwy,bwx] = conditionalmeanmvksdensity(y,x,varargin)
% function [my,bwy,bwx] = conditionalmeanmvksdensity(y,x)
% Computes a conditional mean estimate my for the sample (multivariate) data
% in the n-by-dy matrix y given the sample (multivariate) data in the n-by-d matrix x.
% conditionalmeanmvksdensity evaluates the conditional mean estimate using the multivariate kernel density (mvksdensity) estimation method.
% bwy is the d-element vector containing the bandwidth of the kernel smoothing window in each dimension of y,
% bwx is the d-element vector containing the bandwidth of the kernel smoothing window in each dimension of x, bwy and bwx are the optimal
% bandwidths for estimating normal densities.
% my is the conditional mean value of y given each sample data in x.
% The estimate is based on a ratio of product Gaussian kernel functions, using the window
% parameters (bandwidths bwy and bwx) that are functions of the number of points and dimension in x.
% n is the number of points. d is the number of dimensions of x. 

y = y(:);
[ny,dy] = size(y);
% ny = length(y); % number of data points
% dy = 1; % stochastic dimension for y
[n,d] = size(x);
% n = size(x,1); % number of data points
% d = size(x,2); % stochastic dimension for x
if ~isequal(n,ny)
    error('Wrong number of data points for y and x')
end

bw = @(sig,d,n) sig*(4/((d+2)*n))^(1/(d+4)); % bandwith for kernel density estimation
% sigx = std(x,0); % classical estimate of std
% sigy = std(y,0); % classical estimate of std
sigx = mad(x,1) / 0.6745; % robust estimate of std
sigy = mad(y,1) / 0.6745; % robust estimate of std
bwx = bw(sigx,1+d,n); % multivariate bandwiths
bwy = bw(sigy,1+d,n); % multivariate bandwiths

%% Parfor loop version for computing myx and fx (faster efficiency)
my = zeros(1,n);
if ~verLessThan('matlab','9.2') % introduced in R2017a
    q = parallel.pool.DataQueue;
    afterEach(q, @nUpdateProgressBar);
else
    q = [];
end
textprogressbar('Computing conditional mean: ');
p = 0;
parfor i=1:n
% for i=1:n
    if ~verLessThan('matlab','9.2') % introduced in R2017a
        send(q,i);
    end
    xix = x(i,:);
    myx = meanmvksdensity(y,x,xix,'Bandwidth',bwx);
    if verLessThan('matlab','9.0') % compatibility (<R2016a)
        fx = mymvksdensity(x,xix,'Bandwidth',bwx);
    else
        fx = mvksdensity(x,xix,'Bandwidth',bwx);
    end
    my(:,i) = myx/fx;
end
textprogressbar(' done');

%% For loop version for computing myx (lesser efficiency)
% myx = meanmvksdensity(y,x,x,'Bandwidth',bwx);
% myx = reshape(myx,[1,n]);

%% Parfor loop version for computing myx
% myx = zeros(1,n);
% if ~verLessThan('matlab','9.2') % introduced in R2017a
%     q = parallel.pool.DataQueue;
%     afterEach(q, @nUpdateProgressBar);
% else
%     q = [];
% end
% textprogressbar('Computing joint pdf: ');
% p = 0;
% parfor i=1:n
%     if ~verLessThan('matlab','9.2') % introduced in R2017a
%         send(q,i);
%     end
%     xix = x(i,:);
%     myx(:,i) = meanmvksdensity(y,x,xix,'Bandwidth',bwx);
% end
% textprogressbar(' done');

%% For loop version for computing fx (lesser efficiency)
% if verLessThan('matlab','9.0') % compatibility (<R2016a)
%     fx = mymvksdensity(x,x,'Bandwidth',bwx);
% else
%     fx = mvksdensity(x,x,'Bandwidth',bwx);
% end

% if verLessThan('matlab','9.1') % compatibility (<R2016b)
%     my = bsxfun(@rdivide,myx,fx');
% else
%     my = myx./fx';
% end

function nUpdateProgressBar(~)
p = p+1;
textprogressbar(p/n*100,sprintf('(%d/%d)',p,n));
end

end
