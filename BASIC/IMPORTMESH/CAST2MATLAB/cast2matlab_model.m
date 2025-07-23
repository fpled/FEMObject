function [M,nb]=cast2matlab_model(mode,file,mat)
% function [M,nb]=cast2matlab_model(mode,file,mat)
% MODE : 'UNID', 'PLANE', 'TRID'

fid = fopen(file);
if fid==-1
    error(['Le fichier ' file ' n''existe pas']);
end

[node,elem,nb] = cast2matlab_nodeelem(fid,mode);

if nb.chpo>0
    error(' !!!!!!!!!!   le modele ne doit pas contenir de Champ par point')
end

M = MODEL(mode);
M = addnode(M,node);

if nb.chml>0
    matelem = lirechamp(fid,nb.elem,'MAT');
    matnum = unique(matelem);
else
    matelem = ones(nb.elem,1);
    matnum = 1;
end

for k=1:length(matnum)
    for i=1:length(elem)
        if nargin==3    
            M = addelem(M,getelem(elem{i},find(matelem==matnum(k)),'global'),'mat',getmaterial(mat,matnum(k)));    
        else
            M = addelem(M,getelem(elem{i},find(matelem==matnum(k)),'global'));
        end
    end
end

fclose(fid);

end
