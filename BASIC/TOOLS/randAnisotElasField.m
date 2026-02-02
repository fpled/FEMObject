function C = randAnisotElasField(delta,mL,V)
% function C = randAnisotElasField(delta,mL,V)
% delta : scalar dispersion parameter between 0 and delta_sup =
% sqrt((n3D+1)/(n3D+5)) = sqrt(7/11) = 0.7977 for n3D = 6
% controlling the level of statistical fluctuations of random field G
% mL : n-by-n upper triangular matrix corresponding to the Cholesky factor of n-by-n mean elasticity matrix mC, 
% or n-by-n-by-nx array for n-by-n matrix-valued mean elasticity field mC evaluated at nx points, with n=3 in 2D and n=6 in 3D
% V : n(n+1)/2-by-N-by-nx array containing N samples of centered Gaussian random vector 
% gathering the nU=n(n+1)/2 normalized Gaussian real-valued random fields U_j, j=1,...,n(n+1)/2, evaluated at nx points
% C : n-by-n-by-N-by-nx array containing N samples of n-by-n matrix-valued non-Gaussian random elasticity field C evaluated at nx points

[nU,N,nx] = size(V);
n3D = 6;
switch nU
    case 6
        % dim = 2;
        n = 3;
    case 21
        % dim = 3
        n = 6;
    otherwise
        error('Wrong size for random elasticity field')
end

%% Gaussian stochastic germs for Gaussian (real-valued) and Gamma (positive-valued) random fields
VGauss = triu(ones(n),1);
VGauss = repmat(VGauss,1,1,N*nx);
VGauss(VGauss==1) = V(1:n*(n-1)/2,:,:);

VGamma = eye(n); % VGamma = diag(ones(n,1));
VGamma = repmat(VGamma,1,1,N*nx);
VGamma(VGamma==1) = V(n*(n-1)/2+1:n*(n+1)/2,:,:);
clear V

%% Parameters
sigma = delta/sqrt(n3D+1);
alpha = 1/(2*sigma^2) + (1-(1:n)')/2;
alpha = diag(alpha);
beta = eye(n); % beta = diag(ones(n,1));

%% Parallel (parfor loop) implementation for simulating Gamma (positive-valued) random fields Gam (lesser efficiency)
% Gam = zeros(size(VGamma));
% if ~verLessThan('matlab','9.2') % introduced in R2017a
%     q = parallel.pool.DataQueue;
%     afterEach(q, @nUpdateProgressBar);
% else
%     q = [];
% end
% textprogressbar('Computing Gamma random fields: ');
% p = 0;
% parfor k=1:N*nx
%     if ~verLessThan('matlab','9.2') % introduced in R2017a
%         send(q,k);
%     end
%     Gam(:,:,k) = gaminv(normcdf(VGamma(:,:,k)),alpha,beta);
% end
% textprogressbar(' done');

%% Vectorized implementation for simulating Gamma (positive-valued) random fields Gam (faster efficiency)
alpha = repmat(alpha,1,1,N*nx);
beta = repmat(beta,1,1,N*nx);
Gam = gaminv(normcdf(VGamma),alpha,beta);
clear VGamma

%% Upper triangular matrix-valued random field L and associated transpose Lt
LGauss = sigma * VGauss;
LGamma = sigma * sqrt(2*Gam);
clear VGauss Gam alpha beta
L = LGauss + LGamma;
clear LGauss LGamma
Lt = permute(L,[2,1,3]);

%% Upper triangular matrix or matrix-valued field mL (Cholesky factor of mean elasticity matrix or matrix-valued field mC) and associated transpose mLt 
if ismatrix(mL)
    mLt = mL';
else
    mL = repmat(mL,1,1,1,N);
    mLt = permute(mL,[2,1,3,4]);
end

%% Sequential (for loop) implementation for simulating non-Gaussian matrix-valued random elasticity field C (lesser efficiency)
% C = zeros(n,n,N*nx);
% for k=1:N*nx
%     G = L(:,:,k)'*L(:,:,k);
%     if ismatrix(mL)
%         C(:,:,k) = mL'*G*mL;
%     else
%         C(:,:,k) = mL(:,:,k)'*G*mL(:,:,k);
%     end
% end

%% Parallel (parfor loop) implementation for simulating non-Gaussian matrix-valued random elasticity field C (worst efficiency)
% C = zeros(n,n,N*nx);
% if ~verLessThan('matlab','9.2') % introduced in R2017a
%     q = parallel.pool.DataQueue;
%     afterEach(q, @nUpdateProgressBar);
% else
%     q = [];
% end
% textprogressbar('Computing elasticity random fields: ');
% p = 0;
% parfor k=1:N*nx
%     if ~verLessThan('matlab','9.2') % introduced in R2017a
%         send(q,k);
%     end
%     G = L(:,:,k)'*L(:,:,k);
%     if ismatrix(mL)
%         C(:,:,k) = mL'*G*mL;
%     else
%         C(:,:,k) = mL(:,:,k)'*G*mL(:,:,k);
%     end
% end
% textprogressbar(' done');

%% Vectorized implementation for simulating non-Gaussian matrix-valued random elasticity field C (faster efficiency)
% if ismatrix(mL)
%     mLt = repmat(mLt,1,1,N*nx);
%     mL = repmat(mL,1,1,N*nx);
% end
% if verLessThan('matlab','9.1') % compatibility (<R2016b)
%     G = squeeze(sum(bsxfun(@times,reshape(Lt,[n,n,1,N*nx]),reshape(L,[1,n,n,N*nx])),2));
%     mLtG = squeeze(sum(bsxfun(@times,reshape(mLt,[n,n,1,N*nx]),reshape(G,[1,n,n,N*nx])),2));
%     C = sum(bsxfun(@times,reshape(mLtG,[n,n,1,N*nx]),reshape(mL,[1,n,n,N*nx])),2);
% else
%     G = sum(reshape(Lt,[n,n,1,N*nx]).*reshape(L,[1,n,n,N*nx]),2);
%     mLtG = sum(reshape(mLt,[n,n,1,N*nx]).*reshape(G,[1,n,n,N*nx]),2);
%     C = sum(reshape(mLtG,[n,n,1,N*nx]).*reshape(mL,[1,n,n,N*nx]),2);
% end

if ismatrix(mL)
    mLt = reshape(mLt,[n,n,1]);
    mL = reshape(mL,[1,n,n]);
else
    mLt = reshape(mLt,[n,n,1,N*nx]);
    mL = reshape(mL,[1,n,n,N*nx]);
end
if verLessThan('matlab','9.1') % compatibility (<R2016b)
    G = sum(bsxfun(@times,reshape(Lt,[n,n,1,N*nx]),reshape(L,[1,n,n,N*nx])),2);
    mLtG = sum(bsxfun(@times,mLt,reshape(G,[1,n,n,N*nx])),2);
    C = sum(bsxfun(@times,reshape(mLtG,[n,n,1,N*nx]),mL),2);
else
    G = sum(reshape(Lt,[n,n,1,N*nx]).*reshape(L,[1,n,n,N*nx]),2);
    mLtG = sum(mLt.*reshape(G,[1,n,n,N*nx]),2);
    C = sum(reshape(mLtG,[n,n,1,N*nx]).*mL,2);
end

%% Tianyu Zhang's vectorized implementation for simulating non-Gaussian matrix-valued random elasticity field C
% aa = 1:n*N*nx;
% ai = reshape(aa,n,N*nx);
% bi = repmat(ai,n,1);
% bj = repmat(aa,n,1); % bj = aa(ones(n,1),:);
% ii = reshape(bi,1,numel(bi));
% jj = reshape(bj,1,numel(bj));
% LtMat = sparse(ii,jj,Lt(:));
% 
% LPrm = permute(L,[1,3,2]);
% LMat = reshape(LPrm,n*N*nx,n,1);
% 
% GMat = LtMat * LMat;
% 
% if ismatrix(mL)
%     mLTen = repmat(mL,1,1,N*nx);
% else
%     mLTen = mL;
% end
% mLPrm = permute(mLTen,[1,3,2]);
% mLPrmMat = reshape(mLPrm,n*N*nx,n,1);
% 
% if ismatrix(mL)
%     mLtTen = repmat(mLt,1,1,N*nx);
% else
%     mLtTen = mLt;
% end
% mLtMat = sparse(ii,jj,mLtTen(:));
% 
% BMat = mLtMat * GMat;
% B = permute(reshape(BMat,n,N*nx,n),[1,3,2]);
% BPIMat = sparse(ii,jj,B(:));
% 
% CMat = BPIMat * mLPrmMat;
% C = permute(reshape(CMat,n,N*nx,n),[1,3,2]);

%% Final reshape
C = reshape(C,n,n,N,nx);

%% Text progress bar for parallel implementation
% function nUpdateProgressBar(~)
% p = p+1;
% textprogressbar(p/(N*nx)*100,sprintf('(%d/%d)',p,N*nx));
% end

end
