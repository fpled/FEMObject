function plot_enrichfunmaterial(elem,node,ls,varargin)

connec = calc_conneclocal(elem,node);
issurf = ischarin('surface',varargin);
options = delonlycharin('surface',varargin);
iscourbe = ischarin('courbe',varargin);
options = delonlycharin('courbe',options);

nodecoord = getcoord(node);

if ~isenrich(elem)
   
    if getdim(elem)==1 && iscourbe
    for i=1:size(connec,1)  
    plot(nodecoord(connec(i,:)),zeros(length(connec(i,:)),1));    
    hold on
    end
    else
    patch('Faces',connec,'Vertices',nodecoord,'FaceVertexCData',zeros(size(nodecoord,1),1),options{:});             
    end    
else

    connecenrich = getparam(elem,'connecenrich');
    connecnbddl = getparam(elem,'connecnbddl');
    if getenrichtype(ls)==3
    connecfictitious = getparam(elem,'connecfictitious');
    factnode = ~connecfictitious;
    end

lsval = getvalue(ls);
xnode = getcoord(node,elem);


for e=1:getnbelem(elem);

    connece = connec(e,:);
    nodecoorde = nodecoord(connece,:);
    if isenrich(elem)
    connecenriche = connecenrich(e,:);
    connecnbddle = connecnbddl(e,:);
    if getenrichtype(ls)==3
        factnodee=factnode(e,:)';
    else
        factnodee=ones(3,1);
    end
    end
    
    lse = lsval(connece);
    
    if all(lse<=0)
        subconnecin = 1:getnbnode(elem);
        subconnecout = zeros(0,3);
        xlnodetotal = nodelocalcoord(elem);
        nodecoordeplus = zeros(0,getindim(elem));
        
    elseif all(lse>=0)
        subconnecout = 1:getnbnode(elem);
        subconnecin = zeros(0,3);
        xlnodetotal = nodelocalcoord(elem);
        nodecoordeplus = zeros(0,getindim(elem));
    else
     [subconnecin,subconnecout,xlnodeplus]=lsdivide_oneelem(elem,lse);
     xlnodetotal=[nodelocalcoord(elem);xlnodeplus];    
     nodecoordeplus = getN(elem,xlnodeplus)*nodecoorde; 
    end
    
    nodecoordetotal = [nodecoorde;nodecoordeplus]; 
    
     
    Nuni = getN(elem,xlnodetotal);
    
        
    switch getenrichtype(ls)
        case 1      
        psiin = Nuni*(abs(lse)+lse);
        psiout = Nuni*(abs(lse)-lse);
        case 2
        psiin = -Nuni*lse;
        psiout = Nuni*lse;    
        case 3
            global beta
            if isempty(beta)
                beta=0;
            end
            factnodee=factnode(e,:);
            
            psiin = Nuni*((beta+lse).*factnodee(:));
            psiout = Nuni*((beta-lse).*factnodee(:));
            
    end

if issurf
    nodein = [nodecoordetotal,psiin];
    nodeout = [nodecoordetotal,psiout];       
else
   nodein =  nodecoordetotal;
   nodeout = nodecoordetotal;
end


    if size(subconnecin,1)>0
        
    for k=1:size(subconnecin,1)
        
    if getdim(elem)==1 && iscourbe
    for i=1:size(connec,1)  
    plot(nodein(subconnecin(k,:)),psiin(subconnecin(k,:),:));  
    hold on
    end
    else
    patch('Faces',1:size(subconnecin,2),'Vertices',nodein(subconnecin(k,:),:),'FaceVertexCData',psiin(subconnecin(k,:),:),options{:});             
    end
    end
    
    %patch('Faces',subconnecin,'Vertices',nodein(subconnecin,:),'FaceVertexCData',psiin,options{:});         
    end
    
    if size(subconnecout,1)>0  
        for k=1:size(subconnecout,1)
        if getdim(elem)==1 && iscourbe
    for i=1:size(connec,1)  
    plot(nodeout(subconnecout(k,:)),psiout(subconnecout(k,:),:));    
    hold on
    end
        else
    patch('Faces',1:size(subconnecout,2),'Vertices',nodeout(subconnecout(k,:),:),'FaceVertexCData',psiout(subconnecout(k,:),:),options{:});         
        end
        end

        %patch('Faces',subconnecout,'Vertices',nodeout,'FaceVertexCData',psiout,options{:});
    
    
    end
    
    end
end


end

   