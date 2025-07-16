function [subgaussin,subgaussout] = testsubint(elem,node,ls,ordre)


e = getelem(elem,1);

[elemin,elemout,nodeplus,xnodein,xnodeout ]= lsdivideelem(elem,ls,node);


if isa(ls,'LEVELSET')
lsx = getvalue(ls);
elseif isa(ls,'double')
lsx = ls ;
else
error('')
end


for k=1:length(elemin)
    plot(elemin{k},nodeplus,'FaceColor','y')
end

for k=1:length(elemout)
    plot(elemout{k},nodeplus,'FaceColor','g')
    
end

connece = getconnec(e);

tic
for k=1:1
[subgaussin,subgaussout] = lssubgauss_oneelem(lsx(connece),ordre);
end
toc
plot(POINT(subgaussin.coord),'r*');
plot(POINT(subgaussout.coord),'k*');

sum(subgaussin.w)+sum(subgaussout.w)
axis image