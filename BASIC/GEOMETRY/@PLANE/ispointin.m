function [rep,P] = ispointin(D,P)
% function [rep,P] = ispointin(D,P)

tol = getfemobjectoptions('tolerancepoint');

rep = find(distance(D,P)<tol);

if nargout==2
    P = P(rep);
end
