function [fy,xiy,bwy,bwx] = conditionalksdensity(y,x,varargin)
% function [fy,xiy,bwy,bwx] = conditionalksdensity(y,x)
% function [fy,xiy,bwy,bwx] = conditionalksdensity(y,x,'NumPoints',npts)
% Computes a conditional probability density estimate fy for the sample (univariate) data
% in the vector y given the sample (multivariate) data in the n-by-d matrix x.
% conditionalksdensity evaluates the conditional density estimate at npts equally-spaced points xiy
% covering the range of the data in y. xiy is the set of npts points, npts=100 by default.
% bwy is the bandwidth of the kernel smoothing window for y, bwx is the d-element vector containing the bandwidth of
% the kernel smoothing window in each dimension of x, bwy and bwx are the optimal
% bandwidths for estimating normal densities.
% fy is the npts-by-n matrix of conditional density values of y evaluated at npts points given each sample data in x.
% The estimate is based on a ratio of product Gaussian kernel functions, using the window
% parameters (bandwidths bwy and bwx) that are functions of the number of points and dimension in x.
% n is the number of points. d is the number of dimensions of x. 
% 
% function [fy,xiy,bwy,bwx] = conditionalksdensity(y,x,pts)
% Computes a conditional probability density estimate fy for the sample (univariate) data
% in the vector y, evaluated at the points in the vector pts, given the sample (multivariate) data in the n-by-d matrix x. 
% Here, xiy and pts contain identical values.

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

xiy = [];
if nargin==2 || isempty(varargin)
    m = 100; % number of evaluation points
else
    if ~ischar(varargin{1})
        xiy = varargin{1};
        varargin(1) = [];
        m = size(xiy,1);
    else
        m = getcharin('NumPoints',varargin,100);
    end
end

bw = @(sig,d,n) sig*(4/((d+2)*n))^(1/(d+4)); % bandwith for kernel density estimation
% sigx = std(x,0); % classical estimate of std
% sigy = std(y,0); % classical estimate of std
sigx = mad(x,1) / 0.6745; % robust estimate of std
sigy = mad(y,1) / 0.6745; % robust estimate of std
bwx = bw(sigx,1+d,n); % multivariate bandwiths
bwy = bw(sigy,1+d,n); % multivariate bandwiths
if isempty(xiy)
    foldwidth = 3;
    xiy = linspace(min(y)-foldwidth*bwy,max(y)+foldwidth*bwy,m)';
end

%% Parfor loop version for computing fyx and fx (faster efficiency)
fy = zeros(m,n);
if ~verLessThan('matlab','9.2') % introduced in R2017a
    q = parallel.pool.DataQueue;
    afterEach(q, @nUpdateProgressBar);
else
    q = [];
end
textprogressbar('Computing conditional pdf: ');
p = 0;
parfor i=1:n
    if ~verLessThan('matlab','9.2') % introduced in R2017a
        send(q,i);
    end
    xix = x(i,:);
    if verLessThan('matlab','9.0') % compatibility (<R2016a)
        fyx = mymvksdensity([y,x],[xiy,repmat(xix,m,1)],'Bandwidth',[bwy,bwx]);
    else
        fyx = mvksdensity([y,x],[xiy,repmat(xix,m,1)],'Bandwidth',[bwy,bwx]);
    end
    if verLessThan('matlab','9.0') % compatibility (<R2016a)
        fx = mymvksdensity(x,xix,'Bandwidth',bwx);
    else
        fx = mvksdensity(x,xix,'Bandwidth',bwx);
    end
    fy(:,i) = fyx/fx;
end
textprogressbar(' done');

%% For loop version for computing fyx (lesser efficiency)
% if verLessThan('matlab','9.0') % compatibility (<R2016a)
%     fyx = mymvksdensity([y,x],[repmat(xiy,n,1),kron(x,ones(m,1))],'Bandwidth',[bwy,bwx]);
% else
%     fyx = mvksdensity([y,x],[repmat(xiy,n,1),kron(x,ones(m,1))],'Bandwidth',[bwy,bwx]);
% end
% fyx = reshape(fyx,[m,n]);

%% Parfor loop version for computing fyx
% fyx = zeros(m,n);
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
%     if verLessThan('matlab','9.0') % compatibility (<R2016a)
%         fyx(:,i) = mymvksdensity([y,x],[xiy,repmat(xix,m,1)],'Bandwidth',[bwy,bwx]);
%     else
%         fyx(:,i) = mvksdensity([y,x],[xiy,repmat(xix,m,1)],'Bandwidth',[bwy,bwx]);
%     end
% end
% textprogressbar(' done');

%% For loop version for computing fx (lesser efficiency)
% if verLessThan('matlab','9.0') % compatibility (<R2016a)
%     fx = mymvksdensity(x,x,'Bandwidth',bwx);
% else
%     fx = mvksdensity(x,x,'Bandwidth',bwx);
% end

% if verLessThan('matlab','9.1') % compatibility (<R2016b)
%     fy = bsxfun(@rdivide,fyx,fx');
% else
%     fy = fyx./fx';
% end

function nUpdateProgressBar(~)
p = p+1;
textprogressbar(p/n*100,sprintf('(%d/%d)',p,n));
end

end
