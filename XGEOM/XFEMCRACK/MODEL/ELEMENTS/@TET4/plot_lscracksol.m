function plot_lscracksol(elem,node,q,ls,varargin)
% function plot_lscracksol(elem,node,q,ls,varargin)

if ~isenrich(elem)
    plot_sol(elem,q,varargin{:})
else
    
    connec = calc_conneclocal(elem,node);
    connecenrich = getparam(elem,'connecenrich');
    connecnbddl = getparam(elem,'connecnbddl');
    % conneclsenrichnature = getlsenrichnature(node,connec);
    nodecoord = getcoord(node);
    lssupport=getlssupport(ls);
    lssupportval = getvalue(lssupport);
    xnode = getcoord(node,elem);
    
    issigma =  ischarin('sigma',varargin);
    ampl = getcharin('ampl',varargin);
    varargin=delcharin('ampl',varargin);
    
    if issigma
        ksigma = getcharin('sigma',varargin);
        varargin = delcharin('sigma',varargin);
    end
    
    qelem=localize(elem,q);
    mat = getmaterial(elem);
    
    if issigma
        D = calc_opmat(mat,elem);
        B = calc_B(elem,xnode,[]);
    end
    
    for e=1:getnbelem(elem);
        % eleme= getelem(elem,e);
        connece = connec(e,:);
        nodecoorde = nodecoord(connece,:);
        connecenriche = connecenrich(e,:);
        connecnbddle = connecnbddl(e,:);
        xnodee = xnode(:,:,e);
        
        if issigma
            Be = double(B(:,:,e));
        end
        
        ls1e = lssupportval(connece);
        
        qe = qelem(:,:,e);
        [ae,be] = splitsol(qe,connecnbddle,connecenriche);
        
        if all(ls1e>=0) || all(ls1e<=0)
            if all(ls1e>=0)
                ue = ae + be ;
            elseif all(ls1e<=0)
                ue = ae - be ;
            end
            if issigma
                se = D*(Be*ue(:));% sigma(mat,eleme,xnodee,[1/3,1/3],uet(:));
                se= double(se(ksigma));se = se(:);
                varargin = setcharin('FaceColor',varargin,'flat');
                varargin = setcharin('FaceVertexCData',varargin,se);
            end
            
            patch('Faces',1:3,'Vertices',nodecoorde+ampl*ue',varargin{:});
        else
            % division en sous-elements
            [connecin,connecout,xlnodeplus]=lsdivide_oneelem(ls1e);
            xlnodetotal=[nodelocalcoordtri3();xlnodeplus];
            % nodecoordeplus = calc_x(eleme,nodecoorde,xlnodeplus);
            nodecoordeplus = Ntri3(xlnodeplus)*nodecoorde;
            nodecoordeplus = [nodecoorde;nodecoordeplus];
            Ne = Ntri3(xlnodetotal);
            ue1 = Ne*(ae + be)';
            ue2 = Ne*(ae - be)';
            
            if issigma
                
                ue1t=(ae + be);
                se1 = D*(Be*ue1t(:));
                % se1 = sigma(mat,eleme,xnodee,xlnodetotal,ue1t(:));
                ue2t=(ae - be);
                se2 = D*(Be*ue2t(:));
                % se2 = sigma(mat,eleme,xnodee,xlnodetotal,ue2t(:));
                
                se1=double(se1(ksigma));se1 = repmat(se1(:),size(connecout,1),1);
                se2=double(se2(ksigma));se2 = repmat(se2(:),size(connecin,1),1);
                
                varargin = setcharin('FaceColor',varargin,'flat');
                varargin = setcharin('FaceVertexCData',varargin,se1);
                patch('Faces',connecout,'Vertices',nodecoordeplus+ampl*ue1',varargin{:});
                varargin = setcharin('FaceVertexCData',varargin,se2);
                patch('Faces',connecin,'Vertices',nodecoordeplus+ampl*ue2',varargin{:});
                
            else
                patch('Faces',connecout,'Vertices',nodecoordeplus+ampl*ue1',varargin{:});
                patch('Faces',connecin,'Vertices',nodecoordeplus+ampl*ue2',varargin{:});
                
            end
            
            
        end
        
    end
    
end


function [ae,be] = splitsol(qe,conbddl,coenrich)

ae = zeros(2,3);
be = zeros(2,3);
for i=1:3
    rep = sum(conbddl(1:i-1));
    ae(:,i)=qe(rep+[1:2]);
    if coenrich(i)
        be(:,i)=qe(rep+[3:4]);
    end
end

return


function [ae,be] = splitsolglobal(qe,conbddl,coenrich)

ae = zeros(2,3,size(connbddl,1));
be = zeros(2,3,size(connbddl,1));

for i=1:3
    rep = sum(conbddl(:,1:i-1),2);
    ae(1,i,:)=qe(rep(:,1)+1);
    ae(2,i,:)=qe(rep(:,1)+2);
    repe = find(coenrich(:,i));
    be(1,i,repe)= qe(rep(repe,1)+3);
    be(2,i,repe)= qe(rep(repe,1)+4);
end
return