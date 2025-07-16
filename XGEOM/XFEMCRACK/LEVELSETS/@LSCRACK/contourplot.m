
function varargout=contourplot(c,D,varargin)

if israndom(c)
    warning('la levelset est aleatoire : pas de contourplot possible')
    return
end
contourval=getcharin('value',varargin,0);
c = lseval(c,D);
D = setlevelsets(D,c);
D = lssplitelem(D);


for i=1:length(contourval)
    
    for p=1:getnbgroupelem(D)
        elem = getgroupelem(D,p);
        lstype = getlstype(elem);
        if strcmp(lstype,'cut') || strcmp(lstype,'bicut') || isenrich(elem)
            Dadd = MODEL(getmode(D));
            
            switch lstype
                case {'cut'}
                    [elem,node]=contour(elem,D.node,getvalue(getlssupport(c)),contourval(i));
                    Dadd = addelem(Dadd,elem);
                    Dadd = addnode(Dadd,node);
                    %plot(Dadd,'FaceColor','g','FaceLighting','gouraud','FaceAlpha',.1,'Color','k')
                    %plot(Dadd,'FaceColor','none','FaceLighting','gouraud','FaceAlpha',.1,'Color','k')
                    
                    % keyboard
                case 'bicut'
                    
                    lstip = getlstip(c,getparam(getgroupelem(D,p),'tipnumber'));
                    lssupport = getlssupport(c);
                    [elemcutin,elemcutout,nodeplus,xnodein,xnodeout,ls1value]=...
                        lsdivideelem(elem,lstip,D.node,getvalue(lssupport));
                    
                    %figure
                    %clf
                    %plot(D,'selgroup',3,'FaceColor','none','FaceLighting','gouraud','FaceAlpha',.5,'Color','k')
                    %plot(D,'selgroup',1,'FaceColor','r','FaceLighting','gouraud','FaceAlpha',.5,'Color','k')
                    %plot(D,'selgroup',2,'FaceColor','g','FaceLighting','gouraud','FaceAlpha',.5,'Color','k')
                    %keyboard
                    
                    for k=1:length(elemcutin)
                        [elem,node]=contour(elemcutin{k},nodeplus,ls1value,contourval(i));
                        Dadd = addelem(Dadd,elem);
                        Dadd = addnode(Dadd,node);
                    end
                    
                    %plot(Dadd,'FaceColor','none','FaceLighting','gouraud','FaceAlpha',.5,'Color','b')
                    % Dadd = MODEL(getmode(D));
                    %  for k=1:length(elemcutout)
                    % [elem,node]=contour(elemcutout{k},nodeplus,ls1value,contourval(i));
                    % Dadd = addelem(Dadd,elem);
                    % Dadd = addnode(Dadd,node);
                    % end
                    % plot(Dadd,'FaceColor','w','FaceLighting','gouraud','FaceAlpha',.5,'Color','m')
                    
                    %keyboard
                    
                    
                    
                otherwise
                    %warning('mal programme : prolonge la fissure')
                    %[elem,node]=contour(elem,D.node,getvalue(c.LEVELSETS{1}),contourval(i));
                    %   Dadd = addelem(Dadd,elem);
                    %   Dadd = addnode(Dadd,node);
                    
            end
            
            if length(contourval)==1 | ischarin('Color',varargin)
                col = getcharin('Color',varargin,'r');
                varargin = setcharin('Color',varargin,col);
            else
                col = contourval(i);
            end
            plot(Dadd,'Color',col,varargin{:});
            
        end
    end
end

if nargout>=1
    varargout{1}=node;
    varargout{2}=coseg;
end

%
%
% S.ls = LEVELSETS(c);
% S = lssplitelem(S,'lsenrichtype',2);
% S1 = keepgroupelem(S,'cut');
% contourplot(getlevelset(c,1),S1,varargin{:});
%
% S2 = keepgroupelem(S,'bicut');
%
% if getnbelem(S2)>0
%     warning('fissure mal representee en pointe')
%     contourplot(getlevelset(c,1),S2,varargin{:});
% end
