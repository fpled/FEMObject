function u = ELLIPSOID(varargin)
% function E = ELLIPSOID(E)
% function E = ELLIPSOID(E,c)
% function E = ELLIPSOID(cx,cy,cz,a,b,c)
% function E = ELLIPSOID(cx,cy,cz,a,b,c,vx,vy)
% function E = ELLIPSOID(cx,cy,cz,a,b,c,nx,ny,nz)
% function E = ELLIPSOID(cx,cy,cz,a,b,c,nx,ny,nz,vx,vy)

if nargin==0
    u = ELLIPSOID(0,0,0,1,1,1);
elseif nargin==1
    if isa(varargin{1},'ELLIPSOID')
        u = varargin{1};
    end
elseif nargin==2
    if isa(varargin{1},'ELLIPSE')
        E = varargin{1};
        c = getc(E);
        a = geta(E);
        b = getb(E);
        n = getn(E);
        v = getv(E);
        u.dim = 3;
        u.cx = c(1);
        u.cy = c(2);
        u.cz = c(3);
        u.a  = a;
        u.b  = b;
        u.c  = varargin{2};
        u.nx = n(1);
        u.ny = n(2);
        u.nz = n(3);
        u.vx = v(1);
        u.vy = v(2);
        u.indim = 3;
        
        u = class(u,'ELLIPSOID',GEOMOBJECT(u.dim,u.indim));
    end
elseif nargin==6
    u.dim = 3;
    u.cx = varargin{1};
    u.cy = varargin{2};
    u.cz = varargin{3};
    u.a  = varargin{4};
    u.b  = varargin{5};
    u.c  = varargin{6};
    u.nx = 0;
    u.ny = 0;
    u.nz = 1;
    u.vx = 1;
    u.vy = 0;
    u.indim = 3;
    
    u = class(u,'ELLIPSOID',GEOMOBJECT(u.dim,u.indim));
elseif nargin==8
    u.dim = 3;
    u.cx = varargin{1};
    u.cy = varargin{2};
    u.cz = varargin{3};
    u.a  = varargin{4};
    u.b  = varargin{5};
    u.c  = varargin{6};
    u.nx = 0;
    u.ny = 0;
    u.nz = 1;
    u.vx = varargin{7};
    u.vy = varargin{8};
    u.indim = 3;
    
    u = class(u,'ELLIPSOID',GEOMOBJECT(u.dim,u.indim));
elseif nargin==9
    u.dim = 3;
    u.cx = varargin{1};
    u.cy = varargin{2};
    u.cz = varargin{3};
    u.a  = varargin{4};
    u.b  = varargin{5};
    u.c  = varargin{6};
    u.nx = varargin{7};
    u.ny = varargin{8};
    u.nz = varargin{9};
    u.vx = 1;
    u.vy = 0;
    u.indim = 3;
    
    u = class(u,'ELLIPSOID',GEOMOBJECT(u.dim,u.indim));
elseif nargin==11
    u.dim = 3;
    u.cx = varargin{1};
    u.cy = varargin{2};
    u.cz = varargin{3};
    u.a  = varargin{4};
    u.b  = varargin{5};
    u.c  = varargin{6};
    u.nx = varargin{7};
    u.ny = varargin{8};
    u.nz = varargin{9};
    u.vx = varargin{10};
    u.vy = varargin{11};
    u.indim = 3;
    
    u = class(u,'ELLIPSOID',GEOMOBJECT(u.dim,u.indim));
% else
%     error('Wrong input arguments')
end
end
