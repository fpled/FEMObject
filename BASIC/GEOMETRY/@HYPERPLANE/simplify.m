function D = simplify(D)
% function D = simplify(D)

switch getindim(D)
    case 1
        D = D.P{1};
    case 2
        D = STRAIGHTLINE(D.P{1},D.V{2});
    case 3
        D = PLANE(D.P{1},D.V{1});
end
