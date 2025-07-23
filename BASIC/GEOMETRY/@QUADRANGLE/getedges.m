function [L,N] = getedges(D)
% function [L,N] = getedges(D)

P = getvertices(D);

seg = [1,2;2,3;3,4;4,1];
L = cell(1,4);
for i=1:4
    L{i} = LINE(P{seg(i,1)},P{seg(i,2)});
end

if nargout==2
    N = cell(1,4);
end
