function plot_lsmaterialsol(elem,node,q,ls,varargin)
% function plot_lsmaterialsol(elem,node,q,ls,varargin)

matin = getmaterial(ls);
matout = getmaterial(elem);

switch getlstype(elem)
    case {'in','indomain'}
        elem = setmaterial(elem,matin);
        plot_sol(elem,node,q,varargin{:});
    case {'out'}
        elem = setmaterial(elem,matout);
        plot_sol(elem,node,q,varargin{:});
    case {'cut','touchcut'}
        connec = calc_conneclocal(elem,node);
        
        if isenrich(elem)
            connecenrich = getparam(elem,'connecenrich');
            connecnbddl = getparam(elem,'connecnbddl');
            if getenrichtype(ls)==3
                connecfictitious = getparam(elem,'connecfictitious');
                factnode = ~connecfictitious;
            end
        end
        
        nodecoord = getcoord(node);
        lsval = getvalue(ls);
        xnode = getcoord(node,elem);
        
        issigma =  ischarin('sigma',varargin);
        ampl = getcharin('ampl',varargin);
        varargin=delcharin('ampl',varargin);
        
        if issigma
            ksigma = getcharin('sigma',varargin);
            varargin = delcharin('sigma',varargin);
            varargin = setcharin('FaceColor',varargin,'flat');
        end
        vararginin=varargin;
        vararginout=varargin;
        
        qelem=localize(elem,q);
        
        if issigma
            Din = calc_opmat(matin,elem);
            Dout = calc_opmat(matout,elem);
            [B,detJ,DN] = calc_B(elem,xnode,[]);
        end
        
        for e=1:getnbelem(elem);
            eleme=getelem(elem,e);
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
            
            
            if issigma
                Be = double(B(:,:,e));
                DNe = double(DN(:,:,e));
            end
            
            lse = lsval(connece);
            
            qe = double(qelem(:,:,e));
            
            if isenrich(elem)
                [ae,be] = splitsol(qe,connecnbddle,connecenriche);
            else
                ae = reshape(qe,2,3);
            end
            
            if all(lse<=0) || all(lse>=0)
                if all(lse<=0);choix='in';else;choix='out';
                end
                
                
                if issigma
                    if all(lse<=0);D = Din;else;D=Dout;
                    end
                    if isenrich(elem)
                        Bels= calc_Bls(Be,DNe,[1/3,1/3],lse,choix,getenrichtype(ls),factnodee);
                        se = D*(Bels*[ae(:);be(:)]);
                    else
                        se = D*(Be*qe);
                    end
                    
                    se= double(se(ksigma));se = se(:);
                    varargin = setcharin('FaceVertexCData',varargin,se);
                end
                
                if ~isenrich(elem)
                    ue = ae;
                else
                    Nels = calc_Nls(nodelocalcoordtri3(),lse,choix,getenrichtype(ls),factnodee);
                    ue = Nels*[ae(:);be(:)];
                    ue = reshape(ue,2,3);
                end
                
                patch('Faces',1:3,'Vertices',nodecoorde+ampl*ue',varargin{:});
            else
                % division en sous-elements
                [connecin,connecout,xlnodeplus]=lsdivide_oneelem(lse);
                xlnodetotal=[nodelocalcoordtri3();xlnodeplus];
                % nodecoordeplus = calc_x(eleme,nodecoorde,xlnodeplus);
                nodecoordeplus = Ntri3(xlnodeplus)*nodecoorde;
                nodecoordeplus = [nodecoorde;nodecoordeplus];
                
                if ~isenrich(elem) %%%% non enrichi
                    Ne = Ntri3(xlnodetotal);
                    ue1 = Ne*ae';
                    ue2 = Ne*ae';
                    if issigma
                        ue1t=ae;
                        se1 = Din*(Be*ue1t(:));
                        se1=double(se1(ksigma));
                        se1 = repmat(se1(:),size(connecin,1),1);
                        ue2t=ae;
                        se2 = Dout*(Be*ue2t(:));
                        se2=double(se2(ksigma));
                        se2 = repmat(se2(:),size(connecout,1),1);
                        
                        vararginin = setcharin('FaceVertexCData',varargin,se1);
                        vararginout = setcharin('FaceVertexCData',varargin,se2);
                    end
                    patch('Faces',connecin,'Vertices',nodecoordeplus+ampl*ue1,vararginin{:});
                    patch('Faces',connecout,'Vertices',nodecoordeplus+ampl*ue2,vararginout{:});
                    
                    
                else %%% enrichi
                    
                    for kk=1:size(connecin,1)
                        xlnodek = xlnodetotal(connecin(kk,:),:);
                        xik=sum(xlnodek,1)/3;
                        nodecoordk = nodecoordeplus(connecin(kk,:),:);
                        Nelsk = calc_Nls(xlnodek,lse,'in',getenrichtype(ls),factnodee);
                        uek = Nelsk*[ae(:);be(:)];
                        uek = reshape(uek,2,3);
                        if issigma
                            Belsk = calc_Bls(Be,DNe,xik,lse,'in',getenrichtype(ls),factnodee);
                            se1 = Din*(Belsk*[ae(:);be(:)]);
                            se1=double(se1(ksigma));
                            varargin = setcharin('FaceVertexCData',varargin,se1);
                        end
                        
                        patch('Faces',1:3,'Vertices',nodecoordk+ampl*uek',varargin{:});
                    end
                    
                    for kk=1:size(connecout,1)
                        xlnodek = xlnodetotal(connecout(kk,:),:);
                        xik=sum(xlnodek,1)/3;
                        nodecoordk = nodecoordeplus(connecout(kk,:),:);
                        Nelsk = calc_Nls(xlnodek,lse,'out',getenrichtype(ls),factnodee);
                        uek = Nelsk*[ae(:);be(:)];
                        uek = reshape(uek,2,3);
                        if issigma
                            Belsk = calc_Bls(Be,DNe,xik,lse,'out',getenrichtype(ls),factnodee);
                            se1 = Dout*(Belsk*[ae(:);be(:)]);
                            se1=double(se1(ksigma));
                            varargin = setcharin('FaceVertexCData',varargin,se1);
                        end
                        patch('Faces',1:3,'Vertices',nodecoordk+ampl*uek',varargin{:});
                    end
                    
                    
                    
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



function Bls = calc_Bls(B,DN,xi,ls,choix,enrichtype,factnode)
Nuni = Ntri3(xi);
N = zeros(2,6);
N(1,1:2:end)=Nuni;
N(2,2:2:end)=Nuni;

global beta
switch choix
    case {'out','indomain'}
        switch enrichtype
            case 1
                lsx = Nuni*(abs(ls)-ls);
                dlsx = DN*(abs(ls)-ls);
            case 3
                lsx = Nuni*((beta-ls).*factnode);
                dlsx = DN*((beta-ls).*factnode);
        end
    case 'in'
        switch enrichtype
            case 1
                lsx = Nuni*(abs(ls)+ls);
                dlsx = DN*(abs(ls)+ls);
            case 3
                lsx = Nuni*((beta+ls).*factnode);
                dlsx = DN*((beta+ls).*factnode);
        end
end
Dlsx = zeros(3,2);
Dlsx(1,1) = dlsx(1);
Dlsx(2,2) = dlsx(2);
Dlsx(3,1) = dlsx(2);
Dlsx(3,2) = dlsx(1);
Bls = [B ,lsx*B + Dlsx*N];
return


function Nls = calc_Nls(xi,ls,choix,enrichtype,factnode)
Nuni = Ntri3(xi);
global beta
switch choix
    case {'out','indomain'}
        switch enrichtype
            case 1
                lsx = Nuni*(abs(ls)-ls);
            case 3
                lsx = Nuni*((beta-ls).*factnode);
        end
    case 'in'
        switch enrichtype
            case 1
                lsx = Nuni*(abs(ls)+ls);
            case 3
                lsx = Nuni*((beta+ls).*factnode);
        end
end
Nls = zeros(2*size(Nuni,1),12);

for i=1:size(Nuni,1)
    Nls((i-1)*2+1,1:2:end) = [Nuni(i,:) , lsx(i)*Nuni(i,:)];
    Nls((i-1)*2+2,2:2:end) = [Nuni(i,:) , lsx(i)*Nuni(i,:)];
end

return
