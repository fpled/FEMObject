function u = POINT(varargin)
% u = POINT(coordinates,syscoord)
% syscoord (faculatif) systeme de coordonnees : objet CARTESIAN1D ...

if nargin==0
    u = POINT([0,0]);
elseif nargin==1 && strcmp(class(varargin{1}),'POINT')
    u = varargin{1};
else
    u = struct();
    if isa(varargin{1},'NODE')
        coord = getcoord(varargin{1});
    elseif isa(varargin{1},'VECTOR')
        coord = getcompo(varargin{1})';
    else
        coord = varargin{1};
    end
    
    coord = MYDOUBLEND(coord);
    
    if size(coord,1)>1
        coord = coord';
        coord = reshape(coord,[size(coord,1),1,size(coord,2),sizeND(coord)]);
        coord = coord';
    end
    
    
    if isa(varargin{1},'VECTOR') || isa(varargin{1},'NODE')
        syscoord = getsyscoord(varargin{1});
    else
        if nargin==1
            syscoord = eval(['CARTESIAN' num2str(size(coord,2)) 'D()']);
        else
            syscoord = varargin{2};
        end
    end
    
    u.syscoord = syscoord;
    u = class(u,'POINT',coord,GEOMOBJECT(0,getindim(syscoord)));
    superiorto('MYDOUBLE','MYDOUBLEND','VECTOR','SYSCOORD');
    
end
