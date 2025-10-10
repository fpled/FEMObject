function [rep,P] = ispointin(L,P)
% function [rep,P] = ispointin(L,P)

tol = getfemobjectoptions('tolerancepoint');

v = (L.P{2}-L.P{1})/L.L;
V = P-L.P{1} ;
s = dot(v,V);
rep = find(norm(cross(V,v))<tol & s>=-tol & s<=L.L+tol) ;

if nargout==2
    P = P(rep);
end
