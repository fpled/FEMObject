function [coseg1,node1]=contour(elem,node,nodeval,contourval,varargin)
tolls = getfemobjectoptions('tolerancelevelset');
displayelemtype = 0;

xnode=double(getcoord(node));
dim=size(xnode,2);
[a,connec] = ismember(getconnec(elem),getnumber(node)) ;

dim=size(xnode,2);
con=nodeval(connec)-contourval;
con=reshape(con,size(connec));
repzero=find(abs(con)<tolls);
con(repzero)=0;

seg=zeros(0,2);
for k=1:4
seg=[seg;connec(:,[k,mod(k,4)+1])];
end

segcut=zeros(size(connec,1),4);
for i=1:4
segcut(:,i) = (con(:,i).*con(:,mod(i,4)+1));   
end
%segcut(find(abs(segcut)<tolls))=0;
segcut=sign(segcut);

rep = find(sum(segcut==0 | segcut==-1,2)>1);
segcut = segcut(rep,:);
concut = sign(con(rep,:));

elemin = find(sum(sign(con)==-1 | sign(con)==0,2)==4);

if any(elemin) & displayelemtype
patch('Faces',connec(elemin,:),'Vertices',xnode,'FaceColor','w')
end

% GESTION DES ELEMENTS DE TYPE I
% la level set coupe deux segments 
elem1 = find((sum(segcut==-1,2)==2) );
segcut1 = segcut(elem1,:);
repseg1=find((segcut1==-1)');
repseg1=repseg1(:);
repseg1=reshape(mod(repseg1(:)-1,4)+1,2,length(elem1))';

repseg11=((repseg1(:,1)-1)*size(connec,1)+rep(elem1));
repseg12=mod(repseg11+size(connec,1)-1,numel(connec))+1;
xi1 = con(repseg11)./(con(repseg11)-con(repseg12));
xi1=repmat(xi1,1,dim);
P1 = (1-xi1).*xnode(connec(repseg11),:)+xi1.*xnode(connec(repseg12),:);
repseg11=((repseg1(:,2)-1)*size(connec,1)+rep(elem1));
repseg12=mod(repseg11+size(connec,1)-1,numel(connec))+1;
xi1 = con(repseg11)./(con(repseg11)-con(repseg12));
xi1=repmat(xi1,1,dim);

P2 = (1-xi1).*xnode(connec(repseg11),:)+xi1.*xnode(connec(repseg12),:);


node1 = [P1;P2];
coseg1 = [[1:length(elem1)]',[length(elem1)+1:2*length(elem1)]'];
nbnode=size(node1,1);

if any(elem1) & displayelemtype
patch('Faces',connec(rep(elem1),:),'Vertices',xnode,'FaceColor','b')

end
% GESTION DES ELEMENTS DE TYPE II
% la level set coupe un segment et passe par un point
elem1 = find((sum(segcut==-1,2)==1) & (sum(concut==0,2)==1));
concut1 = concut(elem1,:);
reppo1 = find((concut1==0)');
reppo1 = reppo1(:);
reppo1 = mod(reppo1(:)-1,4)+1;
reppo11 = ((reppo1(:,1)-1)*size(connec,1)+rep(elem1));
P1=xnode(connec(reppo11),:);

segcut1 = segcut(elem1,:);
repseg1=find((segcut1==-1)');
repseg1=repseg1(:);
repseg1=mod(repseg1(:)-1,4)+1;
repseg11=((repseg1(:,1)-1)*size(connec,1)+rep(elem1));
repseg12=mod(repseg11+size(connec,1)-1,numel(connec))+1;
xi1 = con(repseg11)./(con(repseg11)-con(repseg12));
xi1=repmat(xi1,1,dim);
P2 = (1-xi1).*xnode(connec(repseg11),:)+xi1.*xnode(connec(repseg12),:);
if any(elem1)  & displayelemtype
patch('Faces',connec(rep(elem1),:),'Vertices',xnode,'FaceColor','g')


end
node1 = [node1;P1;P2];
coseg1 = [coseg1;...
           nbnode+[[1:length(elem1)]',[length(elem1)+1:2*length(elem1)]']];
nbnode=size(node1,1);

% GESTION DES ELEMENTS DE TYPE III 
% la level passe par 2 noeuds opposes
elem1 = find((sum(segcut==0,2)==4) & (sum(concut==0,2)==2));
concut1 = concut(elem1,:);
reppo1 = find((concut1==0)');

reppo1=reshape(mod(reppo1(:)-1,4)+1,2,length(elem1))';
reppo11 = ((reppo1(:,1)-1)*size(connec,1)+rep(elem1));
reppo12 = ((reppo1(:,2)-1)*size(connec,1)+rep(elem1));
P1=xnode(connec(reppo11),:);
P2=xnode(connec(reppo12),:);

node1 = [node1;P1;P2];
coseg1 = [coseg1;...
           nbnode+[[1:length(elem1)]',[length(elem1)+1:2*length(elem1)]']];
nbnode=size(node1,1);

if any(elem1)  & displayelemtype
patch('Faces',connec(rep(elem1),:),'Vertices',xnode,'FaceColor','y')


end

% GESTION DES ELEMENTS DE TYPE IV
% la level passe par 2 noeuds contigus
elem1 = find((sum(segcut==0,2)==3) & (sum(concut==0,2)==2));
concut1 = concut(elem1,:);
reppo1 = find((concut1==0)');
reppo1=reshape(mod(reppo1(:)-1,4)+1,2,length(elem1))';
reppo11 = ((reppo1(:,1)-1)*size(connec,1)+rep(elem1));
reppo12 = ((reppo1(:,2)-1)*size(connec,1)+rep(elem1));
P1=xnode(connec(reppo11),:);
P2=xnode(connec(reppo12),:);

node1 = [node1;P1;P2];
coseg1 = [coseg1;...
           nbnode+[[1:length(elem1)]',[length(elem1)+1:2*length(elem1)]']];
nbnode=size(node1,1);


if any(elem1)  & displayelemtype
patch('Faces',connec(rep(elem1),:),'Vertices',xnode,'FaceColor','m')

end
% GESTION DES ELEMENTS DE TYPE V
% la level passe par 3 noeuds contigus
elem1 = find((sum(concut==0,2)==3));
concut1 = concut(elem1,:);
reppo1 = find((concut1==0)');
reppo1=reshape(mod(reppo1(:)-1,4)+1,3,length(elem1))';
reppo11 = ((reppo1(:,1)-1)*size(connec,1)+rep(elem1));
reppo12 = ((reppo1(:,2)-1)*size(connec,1)+rep(elem1));
reppo13 = ((reppo1(:,3)-1)*size(connec,1)+rep(elem1));
P1=xnode(connec(reppo11),:);
P2=xnode(connec(reppo12),:);
P3=xnode(connec(reppo13),:);

node1 = [node1;P1;P2;P3];

reps1 = reppo1(:,2)-reppo1(:,1);
reptemp1 = find(reps1==1);
coseg1 = [coseg1;nbnode+[reptemp1,length(elem1)+reptemp1]];
reps1 = reppo1(:,3)-reppo1(:,2);
reptemp1 = find(reps1==1);
coseg1 = [coseg1;nbnode+[length(elem1)+reptemp1,2*length(elem1)+reptemp1]];
reps1 = mod(reppo1(:,1)-reppo1(:,3),4);
reptemp1 = find(reps1==1);
coseg1 = [coseg1;nbnode+[reptemp1,2*length(elem1)+reptemp1]];

nbnode=size(node1,1);


if any(elem1)  & displayelemtype
patch('Faces',connec(rep(elem1),:),'Vertices',xnode,'FaceColor','r')


end
% GESTION DES ELEMENTS DE TYPE V
% la level fait le tour de l'element
elem1 = find((sum(concut==0,2)==4));
for i=1:4
P1 = xnode(connec(rep(elem1),i),:);
P2 = xnode(connec(rep(elem1),mod(i,4)+1),:);
node1 = [node1;P1;P2];
coseg1 = [coseg1;...
           nbnode+[[1:length(elem1)]',[length(elem1)+1:2*length(elem1)]']];
nbnode=size(node1,1);

end    

if any(elem1)  & displayelemtype
patch('Faces',connec(rep(elem1),:),'Vertices',xnode,'FaceColor','w')
end

node1=NODE(node1,1:size(node1,1));
coseg1 = SEG2(node1,1:size(coseg1,1),coseg1);
