function lsplot(M,varargin)

if length(M.ls)==0
    warning('pas de levelset definie pour le modele')
end


display_ = ischarin({'in','out','cut','subcut','cutin','cutout','indomain'},varargin);
if ~any(display_)
    display_([1:3,7])=1;
end
Sin = MODEL(M.mode);
Sin = addnode(Sin,M.node);
Sout=Sin;
Sindomain=Sin;
Scut=Sin;
Scutin=Sin;
Scutout=Sin;

for i=1:M.nbgroupelem
    elem = M.groupelem{i};
    
    switch getlstype(elem)
        case 'indomain'
            Sindomain=addelem(Sindomain,elem,'norenumelem');
        case 'in'
            Sin=addelem(Sin,elem,'norenumelem');
        case 'out'
            Sout=addelem(Sout,elem,'norenumelem');
        case 'cut'
            Scut=addelem(Scut,elem,'norenumelem');
        otherwise
            error(' ')
    end
    
    
    if any(display_(4:6)) & strcmp(getlstype(elem),'cut')
        
        ls = getlevelset(M.ls,getlsnumber(elem));
        [elemcut,repcut] = lsgetelem(elem,ls,'cut',M.node);
        
        if getnbelem(elemcut)>0
            [elemcutin,elemcutout,nodeplus]=lsdivideelem(elemcut,ls,Scutin.node);
            Scutin = addnode(Scutin,nodeplus);
            Scutout = addnode(Scutout,nodeplus);
            for k=1:length(elemcutin)
                Scutin = addelem(Scutin,elemcutin{k},'norenumelem');
            end
            for k=1:length(elemcutout)
                Scutout = addelem(Scutout,elemcutout{k},'norenumelem');
            end
        end
    end
    
end


if display_(1)
    varargin = setcharin('FaceColor',varargin,'y');
    plot(Sin,varargin{:})
end
if display_(2)
    varargin = setcharin('FaceColor',varargin,'g');
    plot(Sout,'out',varargin{:})
end
if display_(3)
    varargin = setcharin('FaceColor',varargin,'m');
    plot(Scut,varargin{:})
end
if display_(4) | display_(5)
    varargin = setcharin('FaceColor',varargin,'y');
    plot(Scutin,varargin{:})
end
if display_(4) | display_(6)
    varargin = setcharin('FaceColor',varargin,'g');
    plot(Scutout,'out',varargin{:})
end
if display_(7)
    varargin = setcharin('FaceColor',varargin,'g');
    plot(Sindomain,varargin{:})
end

if ischarin('enrich',varargin)
    node=M.node;
    for p=1:M.nbgroupelem
        nodep = getnode(node,unique(getconnec(M.groupelem{p})));
        if getlsenrich(M.groupelem{p})>0
            plot(nodep,'Marker','rs')       ;
        else
            if getdim(M.groupelem{p})==1
                plot(nodep,'Marker','k.')       ;
            end
        end
    end
    
end



