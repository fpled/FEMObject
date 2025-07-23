function [segbord,nodebord,xnodeseg,repelembord] = boundary(elem,node)
% function [segbord,nodebord,xnodeseg,repelembord] = boundary(elem,node)

co = getconnec(elem);

seg = zeros(0,2);
for k=1:3
    seg = [seg;co(:,[k,mod(k,3)+1])];
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

numlocal = repmat([1:getnbelem(elem)],1,3);
numlocal(repnu) = [];

numelem = getnumber(elem);
numelem = repmat(numelem(:),3,1);
numelem(repnu) = [];

xnodeseg = nodelocalcoord(elem);
xnodeseg = xnodeseg([1,2,2,3,3,1],:);
xnodeseg = reshape(xnodeseg,[2,3,2]);
xnodeseg = permute(xnodeseg,[1,3,4,2]);
xnodeseg = repmat(xnodeseg,[1,1,getnbelem(elem),1]);
xnodeseg = reshape(xnodeseg,[2,2,getnbelem(elem)*3]);
xnodeseg(:,:,repnu) = [];

segbord = SEG2(nodebord,numelem,segbord,...
    'material',getmaterial(elem),'option','BORD','parent',class(elem));

systri = getsyscoord(elem,numlocal);
sysseg = getsyscoordlocal(segbord);
basetri = getbase(systri);
baseseg = getbase(sysseg);
eY = VECTOR(baseseg(:,1));
if getindim(systri)==3
    eZ = VECTOR(basetri(:,3));
    eX = cross(eY,eZ);
    syscoord = CARTESIAN3D(eX,eY,eZ);
    % syscoord = systri;
else
    eX = rot2D(eY,-pi/2);
    syscoord = CARTESIAN2D(eX,eY);
end

segbord = setsyscoord(segbord,syscoord);

segbord = update_lsdata(segbord,elem);

segbord = setparam(segbord,getparam(elem,numlocal));
