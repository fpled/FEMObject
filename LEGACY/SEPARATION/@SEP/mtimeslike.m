function w = mtimeslike(u,v,fun,TYPE)
% Application de la fonction 'fun' comme un multiplication
% (pas tres clair, hein ?)
% On suppose bien entendu que u et v sont des SEP

[I,J]=ind2sub([u.m,v.m],1:(u.m*v.m));
F=cell(u.m*v.m,u.dim);
F(:,:) = cellfun(fun,u.F(I,:),v.F(J,:),'UniformOutput',0);
alpha = u.alpha(I).*v.alpha(J);
if nargin~=4, SEPfils=str2func(  'SEP' );
else          SEPfils=str2func([ 'SEP' TYPE ]);
end
w=SEPfils(F,alpha);