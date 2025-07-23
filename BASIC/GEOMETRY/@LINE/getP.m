function P = getP(L,k)
% function P = getP(L,k)

if nargin==1
    P = L.P;
else
    P = L.P{k};
end
