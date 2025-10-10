function [rep,P] = ispointin(D,P)
% function [rep,P] = ispointin(D,P)

tol = getfemobjectoptions('tolerancepoint');

c = getcoord(P);

P1 = double(getcoord(D.P1));
P2 = double(getcoord(D.P2));

switch D.dim
    case 1
        rep = find(c(:,1)>=P1(1)-tol & c(:,1)<=P2(1)+tol);
    case 2
        rep = find(c(:,1)>=P1(1)-tol & c(:,1)<=P2(1)+tol & ...
            c(:,2)>=P1(2)-tol & c(:,2)<=P2(2)+tol);
    case 3
        rep = find(c(:,1)>=P1(1)-tol & c(:,1)<=P2(1)+tol & ...
            c(:,2)>=P1(2)-tol & c(:,2)<=P2(2)+tol & ...
            c(:,3)>=P1(3)-tol & c(:,3)<=P2(3)+tol);
end

if nargout==2
    P = P(rep);
end
