function V = getV(u,k)
% function V = getV(u,k)

if nargin==1
   V = u.V;
else
   V = u.V{k};
end
