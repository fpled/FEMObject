function [f,xi,bw] = myksdensity(x,varargin)
% function [f,xi,bw] = myksdensity(x)
% function [f,xi,bw] = myksdensity(x,'Bandwidth',bw,'NumPoints',npts)
% Computes a probability density estimate f for the sample (univariate) data
% in the vector x.
% myksdensity evaluates the density estimate at npts equally-spaced points xi
% covering the range of the data in x. xi is the set of npts points, npts=100 by default.
% bw is the bandwidth of the kernel smoothing window, bw is the optimal
% bandwidth for estimating normal densities. f is the vector of density values.
% The estimate is based on the normal kernel function, using the window
% parameter (bandwidth bw) that is a function of the number of points.
% 
% function [f,xi,bw] = myksdensity(x,pts)
% function [f,xi,bw] = myksdensity(x,pts,'Bandwidth',bw)
% Computes a probability density estimate f for the sample (univariate) data
% in the vector x, evaluated at the points in the vector pts. 
% Here, xi and pts contain identical values.

x = x(:);
[n,d] = size(x);
% n = length(x); % number of data points
% d = 1; % stochastic dimension

xi = [];
if nargin==1 || isempty(varargin)
    m = 100;
    % sig = std(x,0); % classical estimate of std
    sig = mad(x,1) / 0.6745; % robust estimate of std
    bw = sig*(4/((d+2)*n))^(1/(d+4)); % bandwidth
else
    if ~ischar(varargin{1})
        xi = varargin{1};
        varargin(1) = [];
        m = size(xi,1);
    else
        m = getcharin('NumPoints',varargin,100);
    end
    if ischarin('Bandwidth',varargin)
        bw = getcharin('Bandwidth',varargin);
    else
        % sig = std(x,0); % classical estimate of std
        sig = mad(x,1) / 0.6745; % robust estimate of std
        bw = sig*(4/((d+2)*n))^(1/(d+4)); % bandwidth
    end
end

if isempty(xi)
    foldwidth = 3;
    xi = linspace(min(x)-foldwidth*bw,max(x)+foldwidth*bw,m);
end

% mx = mean(x);
% if verLessThan('matlab','9.1') % compatibility (<R2016b)
%     xic = bsxfun(@minus,xi,mx);
%     xc = bsxfun(@minus,x,mx);
% else
%     xc = x - mx;
%     xic = xi - mx;
% end
% % xc = x - repmat(mx,n,1);
% % xic = xi - repmat(mx,1,m);

if verLessThan('matlab','9.1') % compatibility (<R2016b)
    z = bsxfun(@minus,x,xi);
    % z = bsxfun(@minus,xc,xic);
else
    z = x - xi;
    % z = xc - xic;
end
% z = (repmat(x,1,m)-repmat(xi,n,1))/bw;
% z = (repmat(xc,1,m)-repmat(xic,n,1))/bw;
f = normpdf(z);
f = mean(f);
f = f/bw;
