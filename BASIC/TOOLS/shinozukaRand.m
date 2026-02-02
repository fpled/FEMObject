function [V,order,state] = shinozukaRand(x,lcorr,nU,N,varargin)
% function [V,order,state] = shinozukaRand(x,lcorr,nU,N,varargin)
% Generates N independent realizations of nU Gaussian random fields V evaluated 
% at the points in the nx-by-d matrix x with a spatial correlation structure 
% given by the spatial correlation lengths in the d-element vector lcorr 
% using the randomized spectral representation of Shinozuka with the optional 
% name-value pair argument value order for the order (number of terms) 
% of the spectral representation.
% Inputs:
% x is a d-column matrix that specifies the points where the Gaussian 
% random fields are to be evaluated.
% lcorr is a d-element vector specifying the spatial correlation lengths of
% the Gaussian random fields in each spatial dimension or a scalar value for
% all spatial dimensions.
% nx is the number of points x. d is the number of spatial dimensions for x and lcorr.
% nU is the number of Gaussian random fields (1 by default).
% N is the number of independent realizations for each Gaussian random field (1 by default).
% Outputs:
% V is the nx-by-nU-by-N array of double containing the N independent
% realizations of the nU Gaussian random fields evaluatued at the nx points.
% order is the order (number of terms) of the spectral representation.
% state is the structure containing the current random number generator settings after generation.
%
% V = shinozukaRand(x,lcorr,nU,N) returns the nx-by-nU-by-N array of values V
% where the N independent realizations of the nU Gaussian random fields are to be evaluated.
%
% [V,order,state] = shinozukaRand(...) also returns the order 
% (number of terms) of the spectral representation and the structure state 
% containing the current random number generator settings after generation.
%
% [...] = shinozukaRand(...,'PARAM1',val1,'PARAM2',val2,...) specifies parameter 
% name/value pairs to control the order of the spectral representation 
% and the state of the initial random number generator settings before generation. 
% Valid parameters are the following:
%
%  Parameter  Value
%  'order'    A scalar value of the order (number of terms) of the spectral representation. 
%             By default it is computed such that the extent of the domain 
%             contains maximum one half-period in each spatial dimension.
%  'state'    Structure containg the initial random number generator settings
%             before generation.

[nx,d] = size(x);
lcorr = lcorr(:);

if isscalar(lcorr)
    lcorr = repmat(lcorr,d,1);
else
    assert(length(lcorr)==d,'lcorr must be of same length as the number of columns (spatial dimensions) in x')
end

if nargin<3 || isempty(nU)
    nU = 1; % number of Gaussian random fields
end
if nargin<4 || isempty(N)
    N = 1; % number of independent realizations for each Gaussian random field
end
nV = nU*N; % number of independent realizations for all Gaussian random fields

if ischarin('order',varargin)
    order = getcharin('order',varargin); % order of the spectral representation
else
    domainExtent = max(x) - min(x);
    nu = ceil(max(domainExtent./lcorr')); % one-dimensional order nu such that 
    % domainExtent(j) <= period(j)/2 = nu*lcorr(j) for all spatial dimensions j=1,...,d
    nu = 2*floor((nu+1)/2); % ensure one-dimensional order nu is even
    order = nu^d; % d-dimensional order of the spectral representation
end

if ischarin('state',varargin)
    state = getcharin('state',varargin);
    rng(state);
end

X = rand(2+d,order,nV);
Phi = X(1,:,:)*(2*pi); % random phase shifts Phi uniformly distributed on [0,2*pi]
Psi = X(2,:,:); % random variables Psi uniformly distributed on [0,1]
Z = sqrt(-log(Psi)); % random amplitudes Z with values in [0,+Inf[
Tau = (-1+2*X(2+(1:d),:,:)); % random wave numbers Tau uniformly distributed on [-1,1]
clear X Psi

supp = 2*pi./lcorr; % support of power spectral density (PSD) functions

V = zeros(nx,nV);
if ~verLessThan('matlab','9.2') % introduced in R2017a
    qq = parallel.pool.DataQueue;
    afterEach(qq, @nUpdateProgressBar);
else
    qq = [];
end
textprogressbar('Computing Gaussian random fields: ');
p = 0;
parfor i=1:nV
    if ~verLessThan('matlab','9.2') % introduced in R2017a
        send(qq,i);
    end
    tau = Tau(:,:,i);
    q = 1-abs(tau); % triangular function (tri)
    if verLessThan('matlab','9.1') % compatibility (<R2016b)
        k = pi*bsxfun(@rdivide,tau,lcorr); % discretization of wave number k in [-pi/lcorr(j),pi/lcorr(j)] for each spatial dimension j
        s = bsxfun(@times,lcorr,q)/pi; % power spectral density (PSD) functions
    else
        k = pi*tau./lcorr; % discretization of wave number k in [-pi/lcorr(j),pi/lcorr(j)] for each spatial dimension j
        s = lcorr.*q/pi; % power spectral density (PSD) functions
    end
    c = sqrt(prod(supp.*s,1)/order);
    z = Z(:,:,i);
    phi = Phi(:,:,i);
    V(:,i) = shinozukaRand_cpp(c,z,phi,k,x);
end
textprogressbar(' done');

V = reshape(V,nx,nU,N);

if nargout>2
    state = rng;
end

%% Test progress bar
function nUpdateProgressBar(~)
p = p+1;
textprogressbar(p/nV*100,sprintf('(%d/%d)',p,nV));
end

end
