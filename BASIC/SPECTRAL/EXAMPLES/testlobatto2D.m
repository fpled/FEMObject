
n = 10;
m = 6;
if m==1
    La = POLYLAGRANGELOBATTO(n);
else
    La = POLYFELAGRANGE(linspace(-1,1,m),n);    
end
N = getnbpoints(La);
S = SPECTRALMODEL(2,La);
%S = addbc(S,'dirichlet','all');
%
%S = addbc(S,'dirichlet','right');
%S = addbc(S,'dirichlet','left');
%S = addbc(S,'dirichlet','none');

B = calc_matrix(S,1,1);
%
alpha = @(x) exp(-6*(x(:,1).^2+x(:,2).^2));
%alpha = @(x) exp(-20*(x(:,1).^2));
figure(1)
clf
title('loading')
plot(alpha(getpoints(S)),S);

b = calc_vector(S,0,alpha);

u = B\b;

figure(2)
clf
title('solution')
plot(u,S);

%%

%[A,B,C,b]=advec_conv_lobatto(La,alpha,'periodic');

cas=4
%
%
if cas==1
    u0 = @(x) exp(-8*(x(:,1).^2+x(:,2).^2));
    k = 0; 
    vx = -3;
    vy = 0;
elseif cas==2
    k = 0; 
    vx = -3;
    vy = 0;
    u0 = @(x) exp(-20*(x(:,1).^2));
elseif cas==3
    k = .05; 
    v = 15;
    vx = @(x) -x(:,2);
    vy = @(x) x(:,1);
%vx=0;
%vy=1;
    cx = 0.6;
    cy = 0;
    u0 = @(x) exp(-100*((x(:,1)-cx).^2+(x(:,2)-cy).^2));
elseif cas==4
    k = .05; 
    vx = 2;
    vy = 0;
    v=1;
%vx=0;
%vy=1;
    cx = 0;
    cy = 0;
    u0 = @(x) exp(-30*((x(:,1)-cx).^2+(x(:,2)-cy).^2));   
end


u0 = u0(getpoints(S));
u0 = u0(getddlfree(S));
figure(4)
clf
plot(u0,S);

%%
fprintf('calcul A\n')
A = calc_matrix(S,0,0);
fprintf('calcul B\n')
B = calc_matrix(S,1,1);
%C = zeros(size(A));
fprintf('calcul C\n')
C = calc_matrix(S,0,1,vx,vy);
fprintf('fin calcul matrices\n')

%%
%
T = TIMEMODEL(0,.5,150);
N = DGTIMESOLVER(T,0);
ut = dsolve(N,zeros(size(A,1),length(N)),A,k*B+v*C,u0);
zax = [min(min(double(ut))),max(max(double(ut)))];

uref = ut;

ut = setevolparam(ut,'zlim',zax);
ut = setevolparam(ut,'axison',true);
ut = setevolparam(ut,'view',3);

figure(1)

evol(ut,S,'EdgeColor','none')


%%

