function u = VECTEUR(varargin)
% function u = VECTEUR(composantes,syscoord)
% syscoord (faculatif) contient le systeme de coordonnees : 'CARTESIAN' 'CYLINDRIC' 'SPHERIC'
% function u = VECTEUR(POINT)
% function u = VECTEUR(POINT1,POINT2)

if nargin==0
    u = VECTEUR([0;0]);
elseif nargin==1 && isa(varargin{1},'VECTEUR')
    u = varargin{1};
elseif nargin==1 && isa(varargin{1},'SYSCOORD')
    u = VECTEUR(getbase(varargin{1}));
    
elseif nargin==2 && isa(varargin{1},'POINT') && isa(varargin{2},'POINT')
    u = varargin{2}-varargin{1};
else
    u = struct();
    if isa(varargin{1},'POINT')
        compo = getcoord(varargin{1})';
    elseif isa(varargin{1},'double') || isa(varargin{1},'MYDOUBLE')  || isa(varargin{1},'MYDOUBLEND')
        compo = MYDOUBLEND(varargin{1});
        
        if size(compo,2)>1
            compo = reshape(compo,[size(compo,1),1,size(compo,2),sizeND(compo)]);
        end
        
    end
    
    if nargin==1 && isa(varargin{1},'POINT')
        syscoord = getsyscoord(varargin{1});
    elseif nargin==2
        syscoord = varargin{2};
    else
        syscoord = eval(['CARTESIAN' num2str(size(compo,1)) 'D()']);
    end
    
    u.syscoord = syscoord;
    u = class(u,'VECTEUR',compo);
    superiorto('MYDOUBLEND','MYDOUBLE')
end
