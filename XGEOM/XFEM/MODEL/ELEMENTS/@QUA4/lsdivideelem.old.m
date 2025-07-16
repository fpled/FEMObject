function [elemin,elemout,nodeplus,xnodein,xnodeout] = lsdivideelem(elem,ls,node)
LS = getvalue(ls);
connec = getconnec(elem);
connec = [connec,connec(:,1)]';
[a,conneclocal] = ismember(connec,getnumber(node)) ;
lsx = LS(conneclocal);
lsx=reshape(lsx,size(conneclocal));

xnode = getcoord(node,getconnec(elem)');
nodeplus = node;
nbelem = getnbelem(elem);


xnodein = cell(0,1);
xnodeout = cell(0,1);
elemin = cell(0,1);
elemout = cell(0,1);

seg = zeros(2,4,nbelem);
seg(1,:,:) = connec(1:4,:);
seg(2,:,:) = connec(2:5,:);
segcut=sign(lsx(1:4,:).*lsx(2:5,:));

% --------------------------------
% LA LEVELSET PASSE PAR 2 SEGMENTS
% --------------------------------
e1 = find((sum(segcut==-1,1)==2));
if length(e1)>0
connec1 = connec(:,e1);
lsx1 = lsx(:,e1);
[repP1,J] = find(segcut(:,e1)==-1);
repP2 = repP1+1 ;
repP1 = sub2ind(size(lsx1),repP1,J);
repP2 = sub2ind(size(lsx1),repP2,J);
lsxP1 = lsx1(repP1);
lsxP2 = lsx1(repP2);

xi = lsxP1./(lsxP1-lsxP2);
xi=reshape(xi,[2,1,1,length(e1)]);
Pxi = zeros(2,2,4,length(e1));
Pxi(:,1,1,:) = -1+2*xi;
Pxi(:,2,1,:) = -1;
Pxi(:,1,2,:) = 1;
Pxi(:,2,2,:) = -1+2*xi;
Pxi(:,1,3,:) = 1-2*xi;
Pxi(:,2,3,:) = 1;
Pxi(:,1,4,:) = -1;
Pxi(:,2,4,:) = 1-2*xi;

I = find(segcut(:,e1)==-1);
xnodelocal=zerosND(6,2,length(e1)); 
xnodelocal(1:4,:,:) = nodelocalcoord(elem);
xnodelocal(5,:,:)=Pxi(1,:,I(1:2:end));
xnodelocal(5,:,:)=Pxi(1,:,I(1:2:end));
xnodelocal(6,:,:)=Pxi(2,:,I(2:2:end));
xnodelocal(6,:,:)=Pxi(2,:,I(2:2:end));


x=calc_x(elem,xnode(:,:,e1),permute(xnodelocal(5:6,:,:),[4,2,3,1]));
Padd = POINT(x);

numplus = max(getnumber(nodeplus))+[1:numel(Padd)];
nodeplus = addnode(nodeplus,Padd,numplus);

connectri3=zeros(3,0);
connecqua4=zeros(4,0);
xnodetri3 =zerosND(3,2,length(e1));
xnodequa4 =zerosND(4,2,length(e1));

connecplus=[connec(1:4,e1);reshape(numplus,2,length(e1))];

[I,J]=find(segcut(:,e1)==-1);
I=reshape(I,2,length(e1));
rep1 = find(I(1,:)==1 & I(2,:)==2);
rep2 = find(I(1,:)==1 & I(2,:)==3);
rep3 = find(I(1,:)==1 & I(2,:)==4);
rep5 = find(I(1,:)==2 & I(2,:)==3);
rep6 = find(I(1,:)==2 & I(2,:)==4);
rep7 = find(I(1,:)==3 & I(2,:)==4);
connecqua4 = [connecqua4,connecplus([1,5,6,4],rep2)];
connecqua4 = [connecqua4,connecplus([5,2,3,6],rep2)];
connecqua4 = [connecqua4,connecplus([1,2,5,6],rep6)];
connecqua4 = [connecqua4,connecplus([6,5,3,4],rep6)];

jkjhkjhkjh
xnodequa4(:,:,rep1) = xnodelocal([1,4,5,3],:,rep1);
xnodetri3(:,:,rep1) = xnodelocal([4,2,5],:,rep1);
xnodequa4(:,:,rep2) = xnodelocal([4,2,3,5],:,rep2);
xnodetri3(:,:,rep2) = xnodelocal([1,4,5],:,rep2);    
xnodequa4(:,:,rep3) = xnodelocal([1,2,4,5],:,rep3);
xnodetri3(:,:,rep3) = xnodelocal([4,3,5],:,rep3);

num = getnumber(elem,[e1(rep1),e1(rep2),e1(rep3)]);
ssqua=QUA4(nodeplus,num,connecqua4');
sstri=TRI3(nodeplus,num,connectri3');


LS = [LS;zeros(length(numplus),1)];

[ssquain,repquain] = lsgetelem(ssqua,LS,'in',nodeplus);
[ssquaout,repquaout] = lsgetelem(ssqua,LS,'out',nodeplus);

[sstriin,reptriin] = lsgetelem(sstri,LS,'in',nodeplus);
[sstriout,reptriout] = lsgetelem(sstri,LS,'out',nodeplus);

if length(reptriin)>0
elemin=[elemin , {sstriin}];
xnodein=[xnodein , {xnodetri3(:,:,reptriin)}];
end
if length(repquaout)>0
elemout=[elemout , {ssquaout}];
xnodeout=[xnodeout , {xnodequa4(:,:,repquaout)}];
end
if length(repquain)>0
elemin=[elemin , {ssquain}];
xnodein=[xnodein , {xnodequa4(:,:,repquain)}];
end
if length(reptriout)>0
xnodeout=[xnodeout , {xnodetri3(:,:,reptriout)}];
elemout=[elemout , {sstriout}];
end
end


% ----------------------------------------------------------
% LA LEVELSET PASSE PAR UN SOMMET ET COUPE LE SEGMENT OPPOSE
% ----------------------------------------------------------
e2 = find((sum(segcut==-1,1)==1));
if length(e2)>0
connec2 = connec(:,e2);
lsx2 = lsx(:,e2);
[repP1,J] = find(segcut(:,e2)==-1);
repP2 = repP1+1 ;

repP1 = sub2ind(size(lsx2),repP1,J);
repP2 = sub2ind(size(lsx2),repP2,J);
lsxP1 = lsx2(repP1);
lsxP2 = lsx2(repP2);

xi = lsxP1./(lsxP1-lsxP2);
xi=reshape(xi,[1,1,1,length(e2)]);
Pxi = zeros(1,2,3,length(e2));
Pxi(1,1,1,:) = xi;
Pxi(1,1,2,:) = 1-xi;
Pxi(1,2,2,:) = xi;
Pxi(1,2,3,:) = 1-xi;
I = find(segcut(:,e2)==-1);
xnodelocal=zerosND(5,2,length(e2)); 
xnodelocal(1:3,:,:) = nodelocalcoord(elem);
xnodelocal(4,:,:)=Pxi(1,:,I);
xnodelocal(4,:,:)=Pxi(1,:,I);

x=calc_x(elem,xnode(:,:,e2),permute(xnodelocal(5,:,:),[4,2,3,1]));
Padd = POINT(x);

numplus = max(getnumber(nodeplus))+[1:numel(Padd)];
nodeplus = addnode(nodeplus,Padd,numplus);

connectri3=zeros(3,0);
xnodetri3 =zerosND(3,2,length(e2)*2);

connecplus=[connec(1:3,e2);reshape(numplus,1,length(e2))];

[I,J]=find(segcut(:,e2)==-1);
rep1 = find(I==1);
rep2 = find(I==2);
rep3 = find(I==3);
connectri3 = [connectri3,connecplus([1,4,3],rep1),connecplus([4,2,3],rep1)];
connectri3 = [connectri3,connecplus([1,2,4],rep2),connecplus([1,4,3],rep2)];
connectri3 = [connectri3,connecplus([1,2,4],rep3),connecplus([4,2,3],rep3)];

xnodetri3(:,:,rep1) = xnodelocal([1,4,3],:,rep1);
xnodetri3(:,:,rep1+length(e2)) = xnodelocal([4,2,3],:,rep1);
xnodetri3(:,:,rep2) = xnodelocal([1,2,4],:,rep2);
xnodetri3(:,:,rep2+length(e2)) = xnodelocal([1,4,3],:,rep2);
xnodetri3(:,:,rep3) = xnodelocal([1,2,4],:,rep3);
xnodetri3(:,:,rep3+length(e2)) = xnodelocal([4,2,3],:,rep3);

num = getnumber(elem,[e2(rep1),e2(rep2),e2(rep3)]);
sstri=TRI3(nodeplus,[num;num],connectri3');

LS = [LS;zeros(length(numplus),1)];

[sstriin,reptriin] = lsgetelem(sstri,LS,'in',nodeplus);
[sstriout,reptriout] = lsgetelem(sstri,LS,'out',nodeplus);

elemin=[elemin , {sstriin}];
elemout=[elemout , {sstriout}];

xnodein=[xnodein , {xnodetri3(:,:,reptriin)}];
xnodeout=[xnodeout , {xnodetri3(:,:,reptriout)}];

end

%plot(ssqua,nodeplus,'FaceColor','y')
%plot(sstri,nodeplus,'FaceColor','g')


%elemin  = SEG2(nodeplus,getnumber(elem),connecin);
%elemout = SEG2(nodeplus,getnumber(elem),connecout);




if ~isempty(getmaterial(ls))
    matin = getmaterial(ls);
    matout = getmaterial(elem);
else
    matin = getmaterial(elem);
    matout = getmaterial(ls);
end

for k=1:length(elemin)
   elemin{k}=setlstype(elemin{k},'in');
   elemin{k}=setlsnumber(elemin{k},getnumber(ls));
   elemin{k} = setlsenrich(elemin{k},0);
    elemin{k} = setmaterial(elemin{k},matin);

end

for k=1:length(elemout)
if isempty(matout)
   elemout{k}=setlstype(elemout{k},'out');
   elemout{k}=setlsnumber(elemout{k},getnumber(ls));
   elemout{k} = setlsenrich(elemout{k},0);   
end
elemout{k} = setmaterial(elemout{k},getmaterial(elem));
end

