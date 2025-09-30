function P = getvertices(L)
% function P = getvertices(L)

P = L.P;
for i=1:numel(P)
    P{i} = double(getcoord(P{i}));
end