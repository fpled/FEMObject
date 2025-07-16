clear all

f = @(i,j,X,Y) 1/(1+X(i)^2+ sin(2*Y(j))^2);

%% Plot the solution on a small mesh
xplot = linspace(0,1);
yplot = linspace(0,1);

M = numel(xplot);
N = numel(yplot);

Fplot = zeros(M,N);
for i = 1:M
    for j = 1:N
        Fplot(i,j) = f(i,j,xplot,yplot);
    end
end

%% Compute the approximation on a huge mesh
M = 500;
N = 500;
x = linspace(0,1,M);
y = linspace(0,1,N);

clear opts
opts.max_rank = 30;
opts.tol = 1e-12;
opts.err_crit = 'stagnation';

% Pivots obtained with a max/max approach
fprintf('Max/max case\n')
opts.pivot_search = 'max';
[F, ind_i, ind_j] = aca_pp(@(i,j) f(i,j,x,y),M,N,opts);

% Random pivots
fprintf('Random case\n')
opts.pivot_search = 'random_row';
[F_rand, ind_i_rand, ind_j_rand] = aca_pp(@(i,j) f(i,j,x,y),M,N,opts);
%% The error criterion was based on stagnation. The next one is
%% based on entries of the huge matrix
clear opts
opts.max_rank = 15;
opts.tol = 1e-12;
opts.err_crit = 'test_set';
opts.err_card_test_set = min(M,N); % default is M+N;

% Pivots obtained with a max/max approach
fprintf('Max/max case\n')
opts.pivot_search = 'max';
[F, ind_i, ind_j, ind_i_ref, ind_j_ref, F_ref] =...
    aca_pp(@(i,j) f(i,j,x,y),M,N,opts);

% Random pivots
fprintf('Random case\n')
opts.pivot_search = 'random_row';
[F_rand, ind_i_rand, ind_j_rand] = aca_pp(@(i,j) f(i,j,x,y),M,N,opts);


%% Plot
figure(1)
clf
[Xplot, Yplot] = meshgrid(xplot,yplot);
surface(Xplot,Yplot,Fplot,'EdgeColor','none')
colorbar

%%
figure(2)
clf
% Approx on a fine grid
[X,Y] = meshgrid(x,y);
contour(X,Y,F,20)
hold on
plot(x(ind_i_ref),y(ind_j_ref),'ok');
plot(x(ind_i),y(ind_j),'s',...
     'MarkerEdgeColor','r',...
     'MarkerFaceColor','r');
