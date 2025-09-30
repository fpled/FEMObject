function P = getvertices(D)
% function P = getvertices(D)

x1 = double(getcoord(D.P1));
x2 = double(getcoord(D.P2));

switch D.dim
    case 1
        P = cell(1,2);
        P{1} = x1;
        P{2} = x2;
        % P=reshape([P{:}],1,2)';
    case 2
        P = cell(1,4);
        P{1} = [x1(1),x1(2)];
        P{2} = [x2(1),x1(2)];
        P{3} = [x2(1),x2(2)];
        P{4} = [x1(1),x2(2)];
        % P=reshape([P{:}],2,4)';
    case 3
        P = cell(1,8);
        P{1} = [x1(1),x1(2),x1(3)];
        P{2} = [x2(1),x1(2),x1(3)];
        P{3} = [x2(1),x2(2),x1(3)];
        P{4} = [x1(1),x2(2),x1(3)];
        P{5} = [x1(1),x1(2),x2(3)];
        P{6} = [x2(1),x1(2),x2(3)];
        P{7} = [x2(1),x2(2),x2(3)];
        P{8} = [x1(1),x2(2),x2(3)];
        % P=reshape([P{:}],3,8)';
end
