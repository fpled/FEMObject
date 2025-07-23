function P = getvertices(L)
% function P = getvertices(L)

P = L.P;
for i=1:length(P)
    P{i} = double(getcoord(P{i}));
end