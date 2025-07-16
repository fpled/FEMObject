function [elemcontour,nodecontour,contourfield]=contour(elem,node,nodeval,contourval,nodefield)

tolls = getfemobjectoptions('tolerancelevelset');
xnode=double(getcoord(node));
connec = calc_conneclocal(elem,node);

conval=nodeval(connec)-contourval;
conval=reshape(conval,size(connec));
repzero=find(abs(conval)<tolls);
conval(repzero)=0;

repelemcut = find(lsiscut(elem,conval));

seg = [1,2;2,3;1,3;1,4;2,4;3,4];

nodecontour = zeros(0,3);
elemcontour = zeros(0,3);
nbnodes = 0;
nbelemcontour = 0;



for i=1:length(repelemcut)

   e = repelemcut(i);
   connece = connec(e,:);
   xnodee = xnode(connece,:);
   nodeval = conval(e,:);
   segval = nodeval(seg);
   segcut = sign(segval(:,1).*segval(:,2));

   % coupe 3 aretes
   if sum(segcut==-1)==3 
    repseg = find(segcut==-1);
    p=zeros(3,3);
    for k=1:3
    n = seg(repseg(k),:);    
    xi = nodeval(n(1))./(nodeval(n(1))-nodeval(n(2)));  
    p(k,:) = (1-xi)*xnodee(n(1),:)+xi*xnodee(n(2),:);  
    end 
    nodecontour = [nodecontour;p];
    c=1:3;
    elemcontour=[elemcontour;nbnodes+c];

    nbnodes = nbnodes+size(p,1);
    nbelemcontour = nbelemcontour+size(c,1); 
   elseif  sum(segcut==-1)==4
       % coupe 4 aretes
    repseg = find(segcut==-1);
    p=zeros(4,3);
    xi=zeros(1,4);
    for k=1:4
    n = seg(repseg(k),:);       
    xi(k) = nodeval(n(1))./(nodeval(n(1))-nodeval(n(2)));  
    p(k,:) = (1-xi(k))*xnodee(n(1),:)+xi(k)*xnodee(n(2),:);  
    end
    nodecontour = [nodecontour;p];
    c = casefoursegments(repseg);
    %keyboard
    elemcontour=[elemcontour;nbnodes+c];    
    nbnodes = nbnodes+size(p,1);
    nbelemcontour = nbelemcontour+size(c,1); 
    
   elseif sum(nodeval==0)==1 && sum(segcut==-1)==2
       % coupe 2 arete
    repseg = find(segcut==-1);
    repnode = find(nodeval==0);
    p=zeros(3,3);
    p(1,:)=xnodee(repnode,:);
    for k=1:2
    n = seg(repseg(k),:);     
    xi = nodeval(n(1))./(nodeval(n(1))-nodeval(n(2)));  
    p(k+1,:) = (1-xi)*xnodee(n(1),:)+xi*xnodee(n(2),:);  
    end
    c=1:3;
    nodecontour = [nodecontour;p];
    elemcontour=[elemcontour;nbnodes+c];    
    nbnodes = nbnodes+size(p,1);
    nbelemcontour = nbelemcontour+size(c,1); 

    
   elseif sum(nodeval==0)==2 && sum(segcut==-1)==1 
       % coupe 1 arete
    repseg = find(segcut==-1);
    repnode = find(nodeval==0);
    p=zeros(3,3);
    p(1:2,:)=xnodee(repnode,:);
    for k=1
    n = seg(repseg(k),:);     
    xi = nodeval(n(1))./(nodeval(n(1))-nodeval(n(2)));  
    p(k+2,:) = (1-xi)*xnodee(n(1),:)+xi*xnodee(n(2),:);  
    end
     nodecontour = [nodecontour;p];
     c=1:3;
     elemcontour=[elemcontour;nbnodes+c];    
    nbnodes = nbnodes+size(p,1);
    nbelemcontour = nbelemcontour+size(c,1); 

    
   else
       error('pas prevu')
   end
 

%ntemp = NODE(nodecontour(end-size(p,1)+1:end,:),nbnodes-size(p,1)+1:nbnodes);
%elemtemp = TRI3(ntemp,1:size(c,1),elemcontour(end-size(c,1)+1:end,:));

%plot(elemtemp,ntemp,'FaceColor','y','FaceLighting','gouraud','FaceAlpha',.5,'EdgeColor','y')
%keyboard
%plot(elemtemp,ntemp,'FaceColor','none','FaceLighting','gouraud','FaceAlpha',.5,'EdgeColor','k')

end

repelem = find(sum(conval==0,2)>=3);

for i=1:length(repelem)

   e = repelem(i);
   connece = connec(e,:);
   xnodee = xnode(connece,:);
   nodeval = conval(e,:);
   segval = nodeval(seg);
   segcut = sign(segval(:,1).*segval(:,2));

   if sum(nodeval==0)==3 && sum(segcut==-1)==0
    repnode = find(nodeval==0);
    p=zeros(3,3);
    p(1:3,:)=xnodee(repnode,:);  
    nodecontour = [nodecontour;p];
    c=1:3;
    elemcontour=[elemcontour;nbnodes+c];    
    nbnodes = nbnodes+size(p,1);
    nbelemcontour = nbelemcontour+size(c,1); 

   elseif sum(nodeval==0)==4
    repnode = find(nodeval==0);
    p=zeros(4,3);
    p(1:4,:)=xnodee(repnode,:);  
    nodecontour = [nodecontour;p];
    c=[1,2,3;1,2,4;2,3,4;1,3,4];
    elemcontour=[elemcontour;nbnodes+c];    
    nbnodes = nbnodes+size(p,1);
    nbelemcontour = nbelemcontour+size(c,1); 
       warning('cas degenere : fonction nulle sur les 4 sommets du TET4')
   else
       
       error('pas prevu')
   end
   
   
   
end  
    
nodecontour = NODE(nodecontour);
if size(elemcontour,1)>0
elemcontour = TRI3(nodecontour,1:size(elemcontour,1),elemcontour);
else
elemcontour=TRI3();
end

function c=casefoursegments(rep)

[k,cas] = ismember(rep(:)',[1,2,4,6;1,3,5,6;2,3,4,5],'rows');

switch cas
    case 1
        c=[1,2,4;1,4,3];
    case 2
        c=[1,3,2;2,3,4];      
    case 3
        c=[1,2,3;1,3,4];
    otherwise
        error('pas prevu')
end
