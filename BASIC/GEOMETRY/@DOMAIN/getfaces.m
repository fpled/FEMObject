function [F,N] = getfaces(D)
% function [F,N] = getfaces(D)

P = getvertices(D);
switch getdim(D)
    case 3
        faces = [1 2 6 5 ; 2 3 7 6 ; 3 4 8 7 ; 4 1 5 8 ; 1 2 3 4 ; 5 6 7 8];
        F = cell(1,6);
        for i=1:6
            F{i} = PLANE(POINT(P{faces(i,1)}),POINT(P{faces(i,2)}),POINT(P{faces(i,3)}));
        end
        if nargout==2
            N{1} = [0;-1;0];
            N{2} = [1;0;0];
            N{3} = [0;1;0];
            N{4} = [-1;0;0];
            N{5} = [0;0;-1];
            N{6} = [0;0;1];
        end
    case 2
        faces = [1 2 ; 2 3 ; 3 4  ; 4 1 ];
        F = cell(1,4);
        for i=1:4
            F{i} = LINE(POINT(P{faces(i,1)}),POINT(P{faces(i,2)}));
        end
        if nargout==2
            N{1} = [0;-1];
            N{2} = [1;0];
            N{3} = [0;1];
            N{4} = [-1;0];
        end  
    otherwise
        error('pas programme')
        
end