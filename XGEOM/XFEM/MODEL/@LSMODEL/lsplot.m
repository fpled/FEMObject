function lsplot(M,varargin)

if getnblevelsets(M)==0
    warning('pas de levelset definie pour le modele')
end

if ischarin('subcut',varargin)
    if ~ischarin('cutin',varargin)
        varargin = [varargin,'cutin'];
    end
    if ~ischarin('cutout',varargin)
        varargin = [varargin,'cutout'];
    end
end

lstypes = getlstypes();
lstypescolor = getlstypescolor();
display_ = ischarin(lstypes,varargin);
display_ = find(display_);
if ~any(display_)
    lstypes = getlstypes(0);
    lstypescolor = getlstypescolor(0);
else
    lstypes = lstypes(display_);
    lstypescolor = lstypescolor(display_);
end



S = MODEL(getmode(M));
S = addnode(S,getnode(M));

Splot = cell(1,length(lstypes));
Splot(:)={S};

for i=1:getnbgroupelem(M)
    elem = getgroupelem(M,i);
    
    
    if any(ischarin({'cutin','cutout'},varargin)) & (strcmp(getlstype(elem),'cut') || strcmp(getlstype(elem),'bicut'))
        
        ls = getlevelset(M.ls,getlsnumber(elem));
        
        if isa(ls,'LSCRACK')
            [elemcutin,elemcutout,nodeplus]=lscrackdivideelem(elem,ls,getnode(S));
        else
            [elemcutin,elemcutout,nodeplus]=lsdivideelem(elemcut,ls,getnode(S));
        end
        if length(elemcutin)>0
            [temp,jin] = ischarin(getlstype(elemcutin{1}),lstypes);
            Splot{jin}=addnode(Splot{jin},nodeplus);
            for k=1:length(elemcutin)
                Splot{jin} = addelem(Splot{jin},elemcutin{k},'norenumelem');
            end
        end
        if length(elemcutout)>0
            [temp,jout] = ischarin(getlstype(elemcutout{1}),lstypes);
            Splot{jout}=addnode(Splot{jout},nodeplus);
            for k=1:length(elemcutout)
                Splot{jout} = addelem(Splot{jout},elemcutout{k},'norenumelem');
            end
        end
        
    else
        [rep,j] = ischarin(getlstype(elem),lstypes);
        if rep
            Splot{j}=addelem(Splot{j},elem,'norenumelem');
        end
        
    end
    
    
end

leglstypes={};
Handles = [];
for i=1:length(Splot)
    if getnbelem(Splot{i})>0
        varargin = setcharin('FaceColor',varargin,lstypescolor{i});
        Htemp = plot(Splot{i},varargin{:});
        if length(Htemp)>0
            leglstypes = [leglstypes,{lstypes{i}}];
            Handles = [Handles,Htemp(1)];
        end
    end
end


legend(Handles,leglstypes{:});


