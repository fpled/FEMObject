function [L,N] = getedges(D)
% function [L,N] = getedges(D)

P = getvertices(D);
n = numel(P);
seg = [(1:n);circshift(1:n,-1)]';
L = cell(1,n);
for i=1:n-1
    L{i} = LINE(P{seg(i,1)},P{seg(i,2)});
end

if nargout==2
    N = cell(1,n);
end
