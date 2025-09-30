function P = getvertices(D)
% function P = getvertices(D)

P = D.P;
for i=1:numel(P)
    P{i} = double(getcoord(P{i}));
end
