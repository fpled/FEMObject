function plot_lscracksol_tip(elem,node,q,ls,varargin)
% function plot_lscracksol_tip(elem,node,q,ls,varargin)

tipnumber = getparam(elem,'tipnumber') ;
connec = calc_conneclocal(elem,node);
connecenrich = getparam(elem,'connecenrich');
connecnbddl = getparam(elem,'connecnbddl');
conneclsenrichnature = getlsenrichnature(node,connec);
nodecoord = getcoord(node);
lssupport=getlssupport(ls);
lssupportval = getvalue(lssupport);
lstip = getlstip(ls,tipnumber);
lstipval = getvalue(lstip);
xnode = getcoord(node,elem);

issigma =  ischarin('sigma',varargin);
ampl = getcharin('ampl',varargin);
varargin=delcharin('ampl',varargin);

if issigma
    ksigma = getcharin('sigma',varargin);
    varargin = delcharin('sigma',varargin);
    varargin = setcharin('FaceColor',varargin,'flat');
end

qelem=localize(elem,q);
mat = getmaterial(elem);

if issigma
    D = calc_opmat(mat,elem);
    [B,detJ,DN] = calc_B(elem,xnode,[]);
end

% plot(elem,node,'FaceColor','w')
% checke = 12;
% keyboard

for e=1:getnbelem(elem)
    %  if ismember(e,checke)
    %      keyboard
    %  end
    % eleme= getelem(elem,e);
    connece = connec(e,:);
    nodecoorde = nodecoord(connece,:);
    connecenriche = connecenrich(e,:);
    connecnbddle = connecnbddl(e,:);
    conneclsenrichnaturee = conneclsenrichnature(e,:);
    xnodee = xnode(:,:,e);
    
    if issigma
        Be = double(B(:,:,e));
        DNe = double(DN(:,:,e));
    end
    
    ls1e = lssupportval(connece);
    ls2e = lstipval(connece);
    
    qe = qelem(:,:,e);
    [ae,be] = splitsol(qe,connecnbddle,connecenriche);
    be(isnan(be(:)))=0;
    
    if (all(ls1e>=0) || all(ls1e<=0)) && (all(ls2e>=0) || all(ls2e<=0))
        if  all(ls1e>=0)
            Hphi=ones(1,3);
        else
            Hphi=-ones(1,3);
        end
        
        phi2neg = -min(ls2e,0);
        psinode = zeros(1,3);
        
        for i=1:3
            if connecenriche(i) && strncmp(conneclsenrichnaturee{i},'sup',3)
                psinode(i) = Hphi(i);
            elseif connecenriche(i) && strncmp(conneclsenrichnaturee{i},'tip',3)
                psinode(i) = Hphi(i)*phi2neg(i);
            end
        end
        
        
        ue = ae+be.*repmat(psinode,2,1);
        
        
        
        
        if issigma
            se = D*(Be*ue(:));% sigma(mat,eleme,xnodee,[1/3,1/3],uet(:));
            se= double(se(ksigma));se = se(:);
            varargin = setcharin('FaceColor',varargin,'flat');
            varargin = setcharin('FaceVertexCData',varargin,se);
        end
        
        patch('Faces',1:3,'Vertices',nodecoorde+ampl*ue',varargin{:});
        
        
        
    else
        
        % division en sous-elements
        
        [connecin1in2,connecin1out2,connecout1in2,connecout1out2,xlnodeplus]=...
            bilsdivide_oneelem(ls1e,ls2e);
        
        
        subconnec = {connecin1in2,connecin1out2,connecout1in2,connecout1out2};
        xlnodetotal=[nodelocalcoordtri3();xlnodeplus];
        % nodecoordeplus = calc_x(eleme,nodecoorde,xlnodeplus);
        nodecoordeplus = Ntri3(xlnodeplus)*nodecoorde;
        nodecoordeplus = [nodecoorde;nodecoordeplus];
        Ne = Ntri3(xlnodetotal);
        phi2=Ne*ls2e;
        phi2neg = -min(phi2,0);
        
        if issigma
            dphi2 = DNe*ls2e;
        end
        for i=1:4
            psi{i}=zeros(size(Ne,1),3);
        end
        for i=1:3
            if connecenriche(i) && strncmp(conneclsenrichnaturee{i},'sup',3)
                psi{1}(:,i)=-1;
                psi{2}(:,i)=-1;
                psi{3}(:,i)=1;
                psi{4}(:,i)=1;
            elseif connecenriche(i) && strncmp(conneclsenrichnaturee{i},'tip',3)
                psi{1}(:,i)=-1*phi2neg;
                psi{2}(:,i)=-1*phi2neg;
                psi{3}(:,i)=1*phi2neg;
                psi{4}(:,i)=1*phi2neg;
            end
        end
        
        
        de=cell(1,4);
        ue=cell(1,4);
        
        for i=1:4
            % de{i} = (ae+be.*repmat(psinode{i},2,1));
            ue{i} = Ne*ae' + (Ne.*psi{i})*be';
        end
        if issigma
            % warning('ca semble mal fait')
        end
        for i=1:4
            if size(subconnec{i},1)>0
                if issigma
                    se=zeros(size(subconnec{i},1),1);
                    for k=1:size(subconnec{i},1)
                        localcenter = xlnodetotal(subconnec{i}(k,:),:);
                        localcenter = sum(localcenter,1)/3;
                        phi1center = Ntri3(localcenter)*ls1e;
                        phi2center = Ntri3(localcenter)*ls2e;
                        Hphicenter = (phi1center>=0) - (phi1center<0);
                        dphi2center = DNe*ls2e;
                        
                        Ne = zeros(2,6);
                        Ne(1,1:2:end)= Ntri3(localcenter);
                        Ne(2,2:2:end)= Ntri3(localcenter);
                        
                        Bels = zeros(size(Be));
                        for ii=1:3
                            repddli = 2*(ii-1)+[1:2];
                            Bi = Be(:,repddli);
                            Ni = Ne(:,repddli);
                            if connecenriche(ii) && strncmp(conneclsenrichnaturee{ii},'sup',3)
                                psicenter = Hphicenter;
                                Bels(:,repddli)=psicenter*Bi;
                            elseif connecenriche(ii) && strncmp(conneclsenrichnaturee{ii},'tip',3)
                                psicenter = Hphicenter.*(-min(phi2center,0));
                                dpsicenter = Hphicenter*(-dphi2center).*(phi2center<=0);
                                Dpsicenter = zeros(3,2);
                                Dpsicenter(1,1) = dpsicenter(1);
                                Dpsicenter(2,2) = dpsicenter(2);
                                Dpsicenter(3,1) = dpsicenter(2);
                                Dpsicenter(3,2) = dpsicenter(1);
                                Bels(:,repddli)=psicenter*Bi+Dpsicenter*Ni;
                            end
                        end
                        setemp = D*(Be*ae(:)+Bels*be(:));
                        se(k)=double(setemp(ksigma));
                    end
                    varargin = setcharin('FaceVertexCData',varargin,se);
                end
                
                patch('Faces',subconnec{i},'Vertices',nodecoordeplus+ampl*ue{i},varargin{:});
                
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
