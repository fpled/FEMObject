function [segbord,nodebord,xnodeseg,repelembord] = boundary(elem,node)
% function [segbord,nodebord,xnodeseg,repelembord] = boundary(elem,node)

co = getconnec(elem);

seg = zeros(0,2);
for k=1:4
    seg = [seg;co(:,[k,mod(k,4)+1])];
end
nseg = size(seg,1);
segsort = sort(seg,2);
[segu,repu,rep2] = unique(segsort,'rows');
repnu = setdiff([1:nseg],repu);
segnu = segsort(repnu,:);
repnu = find(ismember(segsort,segnu,'rows'));

segbord = seg;
segbord(repnu,:) = [];

nodebord = getnode(node,(unique(segbord)));

numlocal = repmat([1:getnbelem(elem)],1,4);
numlocal(repnu) = [];

numelem = getnumber(elem);
numelem = repmat(numelem(:),4,1);
numelem(repnu) = [];

xnodeseg = nodelocalcoord(elem);
xnodeseg = xnodeseg([1,2,2,3,3,4,4,1],:);
xnodeseg = reshape(xnodeseg,[2,4,2]);
xnodeseg = permute(xnodeseg,[1,3,4,2]);
xnodeseg = repmat(xnodeseg,[1,1,getnbelem(elem),1]);
xnodeseg = reshape(xnodeseg,[2,2,getnbelem(elem)*4]);
xnodeseg(:,:,repnu) = [];

segbord = SEG2(nodebord,(numelem),segbord,'material',getmaterial(elem),'option','BORD','parent',class(elem));

sysqua = getsyscoord(elem,numlocal);
sysseg = getsyscoordlocal(segbord);
basequa = getbase(sysqua);
baseseg = getbase(sysseg);
eY = VECTOR(baseseg(:,1));

if getindim(sysqua)==3
    eZ = VECTOR(basequa(:,3));
    eX = cross(eY,eZ);
    syscoord = CARTESIAN3D(eX,eY,eZ);
    syscoord = sysqua;
else
    eX = rot2D(eY,-pi/2);
    syscoord = CARTESIAN2D(eX,eY);
end

segbord = setsyscoord(segbord,syscoord);

segbord = update_lsdata(segbord,elem);
