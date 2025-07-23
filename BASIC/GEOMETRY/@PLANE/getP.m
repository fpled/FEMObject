function P = getP(u,k)
% function P = getP(u,k)

if nargin==1
    P = u.P;
else
    P = u.P{k};
end
