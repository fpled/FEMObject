function [node,elem,nb] = gmsh2femobject_nodeelem(fid,indim)
% function [node,elem,nb] = gmsh2femobject_nodeelem(fid,indim)


while ~strcmp(fgetl(fid),'$Nodes')
end

nbnode = fscanf(fid,'%f',1);
node = fscanf(fid,'%f %f %f %f',[4,nbnode])';

while ~strcmp(fgetl(fid),'$Elements')
end
nbelem = fscanf(fid,'%f',1);

nbgroup = 0;
elem = cell(1,0);
elemtypetemp = '';

repelemingroup = zeros(nbelem,2);

% tableelem = fscanf(fid,'%f',[33*nbelem])';
tableelem = fscanf(fid,'%f')';
physicaltemp = 0;
geometricaltemp = 0;
partition = 0;

c = 0;
for p=1:nbelem
    elemnum = tableelem(c+1);
    c = c+1;
    [elemtype,nbnodep,dimp] = gmshelemtype(tableelem(c+1));
    c = c+1;
    nbtags = tableelem(c+1);
    c = c+1;
    tags = tableelem(c+(1:nbtags));
    physical = tags(1);
    geometrical = tags(2);
    if nbtags>=3
        partition = tags(3);
    else
        partition = geometrical;
    end
    c = c+nbtags;
    connec = tableelem(c+(1:nbnodep));
    c = c+nbnodep;
    
    if ~strcmp(elemtype,elemtypetemp) || (physical~=physicaltemp) || (geometrical~=geometricaltemp) || (partitiontemp~=partition)
        nbelemingroup = 1;
        nbgroup = nbgroup+1;
        elem{nbgroup}.connec = zeros(nbelem,nbnodep);
        elem{nbgroup}.connec(nbelemingroup,:) = connec ;
        elem{nbgroup}.number = zeros(nbelem,1);
        elem{nbgroup}.number(nbelemingroup,:) = elemnum;
        elem{nbgroup}.dim = dimp;
        elem{nbgroup}.type = elemtype ;
        elem{nbgroup}.physicalentity = physical ;
        elem{nbgroup}.geometricalentity = geometrical ;
        elem{nbgroup}.partition = partition ;
        elemtypetemp = elemtype;
        physicaltemp = physical;
        geometricaltemp = geometrical;
        partitiontemp = partition ;
    else
        nbelemingroup = nbelemingroup+1;
        elem{nbgroup}.connec(nbelemingroup,:) = connec;
        elem{nbgroup}.number(nbelemingroup) = elemnum;
    end
    elem{nbgroup}.nbelemingroup = nbelemingroup;
    
    repelemingroup(p,:) = [nbgroup,nbelemingroup];
end

nb.node = nbnode;
nb.elem = nbelem;
nb.repelemingroup = repelemingroup;

switch indim
    case 1
        node = NODE(node(:,2),node(:,1));
    case 2
        node = NODE(node(:,2:3),node(:,1));
    case 3
        node = NODE(node(:,2:4),node(:,1));
    otherwise
        error('')
end

for i=1:length(elem)
    elem{i}.connec = elem{i}.connec(1:elem{i}.nbelemingroup,:);
    elem{i}.number = elem{i}.number(1:elem{i}.nbelemingroup,:);
    
    % if elem{i}.dim>=1
    physicalentity =    elem{i}.physicalentity;
    geometricalentity = elem{i}.geometricalentity;
    partition = elem{i}.partition;
    elem{i} = eval([elem{i}.type '(node,elem{i}.number,elem{i}.connec)']);
    
    elem{i} = setparam(elem{i},'physicalentity',physicalentity);
    elem{i} = setparam(elem{i},'geometricalentity',geometricalentity);
    elem{i} = setparam(elem{i},'partition',partition);
    %  end
end


function [name,nbnode,dim] = gmshelemtype(number)
% function [name,nbnode,dim] = gmshelemtype(number)

switch number
    case 1
        nbnode = 2;
        name = 'SEG2';
        dim = 1;
    case 2
        nbnode = 3;
        name = 'TRI3';
        dim = 2;
    case 3
        nbnode = 4;
        name = 'QUA4';
        dim = 2;
    case 4
        nbnode = 4;
        name = 'TET4';
        dim = 3;
    case 5
        nbnode = 8;
        name = 'CUB8'; % 'HEXA8'
        dim = 3;
    case 6
        nbnode = 6;
        name = 'PRI6'; % 'WEDGE6'
        dim = 3;
    case 7
        nbnode = 5;
        name = 'PYR5';
        dim = 3;
    case 8
        nbnode = 3;
        name = 'SEG3';
        dim = 1;
    case 9
        nbnode = 6;
        name = 'TRI6';
        dim = 2;
    case 10
        nbnode = 9;
        name = 'QUA9';
        dim = 2;
    case 11
        nbnode = 10;
        name = 'TET10';
        dim = 3;
    case 12
        nbnode = 27;
        name = 'CUB27'; % 'HEXA27'
        dim = 3;
    case 13
        nbnode = 18;
        name = 'PRI18'; % 'WEDGE18'
        dim = 3;
    case 14
        nbnode = 14;
        name = 'PYR14';
        dim = 3;
    case 15
        nbnode = 1;
        name = 'ELEMPOINT';
        dim = 0;
    case 16
        nbnode = 8;
        name = 'QUA8';
        dim = 2;
    case 17
        nbnode = 20;
        name = 'CUB20'; % 'HEXA20'
        dim = 3;
    case 18
        nbnode = 15;
        name = 'PRI15'; % 'WEDGE15'
        dim = 3;
    case 19
        nbnode = 13;
        name = 'PYR13';
        dim = 3;
    otherwise
        error('pas defini')
end

return
