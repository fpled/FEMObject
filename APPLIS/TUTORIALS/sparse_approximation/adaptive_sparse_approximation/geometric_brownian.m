function X = geometric_brownian(x,c,s,T,n)
% function X = geometric_brownian(x,c,s,T,n)
% x: N×d matrix of coefficients
% c, s: scalars
% T, n: final time and number of time steps
% X: N×(n+1) matrix of samples
% Construct a Brownian‐like path B via a truncated sine expansion,
% then run an Euler discretization of dX = c X dt + s X dB.

t  = linspace(0,T,n+1); % time discretization vector t of size 1x(n+1)
dt = T/n;               % time increment dt
[N,d] = size(x);        % number of samples N, parametric dimension d

% Build matrix S of size d×(n+1) so that S(i,k) = sqrt(2)/pi * sin(pi*freq(i)*t(k))/freq(i) for i=1:d, k=1:n+1
freq  = ((1:d)-0.5)';     % frequency vector of size dx1 so that freq(i) = i-0.5 for i=1:d
S = sin(pi*freq*t)./freq; % sine matrix S of size dx(n+1)
S = (sqrt(2)/pi)*S;       % scale each row by sqrt(2)/pi    

% Build truncated sine expansion matrix (Brownian‐like path) B of size N×(n+1) as B = x*S
B = x*S; % Brownian‐like path B of size N×(n+1)

% Compute increment matrix dB of N×n so that dB(:,k) = B(:,k+1)-B(:,k) for k=1:n
dB = diff(B,1,2); % increment matrix dB of size N×n
% dB = B(:,2:end) - B(:,1:end-1);

% Pre‐compute Euler increment matrix Inc so that Inc(:,k) = 1 + c*dt + s*dB(:,k) for k=1:n
Inc = 1 + c*dt + s*dB; % Euler increment matrix Inc of size N×n

% Compute sample matrix X of size N×(n+1) so that X(:,1)=1 and X(:,k+1) = X(:,k) .* Inc(:,k) for k=1:n
% by starting at 1 and taking cumulative products
X = [ones(N,1),cumprod(Inc,2)]; % sample matrix X of size N×(n+1)

end
