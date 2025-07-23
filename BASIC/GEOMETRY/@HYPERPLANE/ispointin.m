function [rep,P] = ispointin(D,P)
% function [rep,P] = ispointin(D,P)

rep = ispointin(simplify(D),P);

if nargout==2
    P = P(rep);
end
