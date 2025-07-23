function [node,elem,nb]=cast2matlab_nodeelem(fid,mode)

A=fscanf(fid,'%f %f %f %f %f',[1,5]);
nbnode=A(1);
nbelem=A(2);
nbchpo=A(3);
nbchml=A(4);

node=zeros(nbnode,4);
nbgroup = 0;
elem=cell(1,0);
node(:,[1,2:4])=fscanf(fid,'%f %f %f %f',[4,nbnode])';
typeelemptemp = '';
matnum = 0;
repelemingroup=zeros(nbelem,2);
for p=1:nbelem

    temp1=fscanf(fid,'%f %f',[2,1])';
    typeelemp=fscanf(fid,'%s',1);
    switch typeelemp
    case 'line'
        elemtype='SEG2';
        connec=fscanf(fid,'%f %f ',[2,1])';
    case 'tri'
        elemtype='TRI3';
        connec=fscanf(fid,'%f %f %f',[3,1])';
    case 'quad'
        elemtype='QUA4';
        connec=fscanf(fid,'%f %f %f %f',[4,1])';

    otherwise
        error('type element non reconnu');        
    end
    if ~strcmp(typeelemp,typeelemptemp)
        nbgroup = nbgroup+1;
        elem{nbgroup}.connec = connec ; 
        elem{nbgroup}.number = p;
    else
        elem{nbgroup}.connec=[elem{nbgroup}.connec;connec];
        elem{nbgroup}.number = [elem{nbgroup}.number;p];    
    end
    elem{nbgroup}.type=elemtype ;
    typeelemptemp=typeelemp;
    repelemingroup(p,:)=[nbgroup,size(elem{nbgroup}.connec,1)];
end

for i=1:length(elem)
    elem{i}.nbelem = size(elem{i}.connec,1); 
end

nb.node=nbnode;
nb.elem=nbelem;
nb.chpo=nbchpo;
nb.chml=nbchml;
nb.repelemingroup=repelemingroup;

switch mode
case 'UNID'
    node = NODE(node(:,2),node(:,1));
case 'PLANE'
    node = NODE(node(:,2:3),node(:,1));        
case 'TRID'
    node = NODE(node(:,2:4),node(:,1));
otherwise
    error('')
end
for i=1:length(elem)
    elem{i}=eval([elem{i}.type '(node,elem{i}.number,elem{i}.connec)']);   
end
