function [u,v] = samesize(u,v)
% function [u,v] = samesize(u,v)

n1 = size(u,2);
n2 = size(v,2);
if n1==1
    u = repmat(u,1,n2);
end
if n2==1
    v = repmat(v,1,n1);
end
if n1>1 && n2>1 && n1~=n2
    error('les vecteurs doivent etre de la meme taille')
end


