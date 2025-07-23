function [P,repelim,repreplace] = unique(P)
% function [P,rep1,rep2] = unique(P)
% elimine les points de P distants de moins que la tolerance eps
% P(rep1) sont les noeuds apres elimination
% P(rep2) sont les noeuds supprimes
%
% See also POINT/distance, POINT/minus, POINT/mtimes, POINT/mrdivide, POINT/ne, POINT/eq,
% POINT/plus, POINT/uminus, POINT/norm, VECTOR/ne

prec = getfemobjectoptions('tolerancepoint');
coord = permute(double(getcoord(P)),[3 2 1]);
[coordu,a,b] = unique(funprec(coord,prec),'rows');
repelim = setdiff(1:size(coord,1),a);
repreplace = a(b(repelim));
coord(repelim,:)=[];

P = POINT(permute(coord,[3,2,1]));


function coord = funprec(coord,prec)

coord = round(coord/prec)*prec;

return
