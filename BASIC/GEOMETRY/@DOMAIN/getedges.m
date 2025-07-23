function [L,N] = getedges(D)
% function [L,N] = getedges(D)

P = getvertices(D);
switch getdim(D)
    case 2
        seg = [1,2;2,3;3,4;4,1];
        L = cell(1,4);
        for i=1:4
            L{i} = LINE(P{seg(i,1)},P{seg(i,2)});
        end
        if nargout==2
            N{1} = [0;-1];
            N{2} = [1;0];
            N{3} = [0;1];
            N{4} = [-1;0];
        end
    case 3
        seg = [1,2;2,3;3,4;4,1];
        seg = [seg;seg+4];
        seg = [seg;[1,5;2,6;3,7;4,8]];
        L = cell(1,12);
        for i=1:12
            L{i} = LINE(P{seg(i,1)},P{seg(i,2)});
        end
        if nargout==2
            error('normale pas definie sur les edges')
        end
    otherwise
        error('pas programme')
        
end