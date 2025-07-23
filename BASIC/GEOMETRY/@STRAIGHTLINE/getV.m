function V = getV(D,k)
% function V = getV(D,k)

if nargin==1
   V = D.V;
else
   V = D.V{k};
end
