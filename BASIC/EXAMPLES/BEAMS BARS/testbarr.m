function testbarr(n)
% function testbarr(n)
% n : number of elements (10 by default)

if nargin<1 || isempty(n)
    n = 10;
end

L = 1;
d = 0.1;
% E = inline('2-x.^2','x');
E = @(x) 2-x.^2;
MAT = ELAS_ISOT('E',E,'NU',0.3,'RHO',1,'S',1);

P1 = POINT([0]);
P2 = POINT([L]);
S1 = mesh(LINE(P1,P2),n,MAT);
S = MODEL('UNID');
S = addelem(S,'BARR',S1,'mat',MAT);

S = final(S);
S = addcl(S,P1,'U',0);

K = calc_rigi(S,'intorder',1);

ampl = 0.1; 

f = nodalload(S,P2,'FX',1);    

q = K\f;
q = unfreevector(S,q);

figure(2)
clf
x = getcoord(S.node);
plot(x,double(q),'-','Marker','.')
hold on
%solex = log(1+x); % solution pour E(x)=1+x
%solex = atan(x);% solution pour E(x)=1+x^2
solex = log((sqrt(2)+x)./(sqrt(2)-x))/2/sqrt(2);% solution pour E(x)=2-x^2
dsolex = 1/2/sqrt(2)*(1./(x+sqrt(2))+1./(sqrt(2)-x));
dq = calc_epsilon(S,q,'smooth');
dq = double(dq);
% dq = double(dq{1});
% solex = x;
figure(1)
hold off
plot(x,solex,'r-','LineWidth',2)
hold on
plot(x,double(q),'b-')
figure(2)
hold off
plot(x,dsolex,'r-','LineWidth',2)
hold on
plot(x,dq(:),'b-')

error1 = norm(solex(end)-double(q(end)))/norm(solex(end))
error2 = norm(dsolex(end)-double(dq(end)))/norm(dsolex(end))

% E = RANDFIELD(RFMARGINAL(RVUNIFORM(0.5,1.5)),EXPCORREL(1));
% E = KL(E,10,S);
