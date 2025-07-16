function varargout=contourplot(ls,D,varargin)
%function varargout=contourplot(ls,D,varargin)
% D: MODEL
%function varargout=contourplot(ls,D,'value',s,varargin)
% s : isovaleur (0 par defaut)

if israndom(ls)
    warning('la levelset est aleatoire : pas de contourplot possible')
    return
end

ls = actualise(ls,D,varargin{:});

if length(ls)>1
    error('ls est une multilevelset')
end

contourval=getcharin('value',varargin,0);
ls=lseval(ls,D);

options = patchoptions(getdim(D)-1,varargin{:});

for i=1:length(contourval)
    for p=1:D.nbgroupelem
        Dadd = MODEL(getmode(D));
        
        [elem,node]=contour(getgroupelem(D,p),getnode(D),getvalue(ls),contourval(i));
        
        if isa(elem,'ELEMENTGEOM')
            dimelem = getdim(elem);
            if dimelem>0
                Dadd = addelem(Dadd,elem);
                Dadd = addnode(Dadd,node);
            else
                error('pas prevu')
            end
        else
            dimelem=0;
            Dadd=node;
        end
        if dimelem>0
            if length(contourval)==1 || ischarin('Color',varargin)
                col = getcharin('Color',varargin,getcharin('EdgeColor',varargin,'r'));
                varargin = setcharin('Color',varargin,col);
            else
                col = contourval(i);
                varargin = setcharin('Color',varargin,col);
            end
            plot(Dadd,varargin{:});
        else
            plot(POINT(node),varargin{:});
        end
        
    end
    
    
end



if nargout>=1
    varargout{1}=node;
    varargout{2}=coseg;
end