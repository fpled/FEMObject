function P = getP(D,k)
% function P = getP(D,k)

if nargin==1
    P = D.P;
else
    P = D.P{k};
end
