function u = ELLIPSE(varargin)
% function E = ELLIPSE(C)
% function E = ELLIPSE(E)
% function E = ELLIPSE(E,a,b)
% function E = ELLIPSE(cx,cy,a,b)
% function E = ELLIPSE(cx,cy,a,b,vx,vy)
% function E = ELLIPSE(cx,cy,cz,a,b)
% function E = ELLIPSE(cx,cy,cz,a,b,vx,vy)
% function E = ELLIPSE(cx,cy,cz,a,b,nx,ny,nz)
% function E = ELLIPSE(cx,cy,cz,a,b,nx,ny,nz,vx,vy)

if nargin==0
    u = ELLIPSE(0,0,1,1);
elseif nargin==1
    if isa(varargin{1},'ELLIPSE')
        u = varargin{1};
    elseif isa(varargin{1},'CIRCLE')
        C = varargin{1};
        c = getc(C);
        n = getn(C);
        v = getv(C);
        r = getr(C);
        u.dim = getdim(C);
        u.cx = c(1);
        u.cy = c(2);
        u.cz = c(3);
        u.a  = r;
        u.b  = r;
        u.nx = n(1);
        u.ny = n(2);
        u.nz = n(3);
        u.vx = v(1);
        u.vy = v(2);
        u.indim = getindim(C);
        
        u = class(u,'ELLIPSE',GEOMOBJECT(u.dim,u.indim));
    end
elseif nargin==3
    if isa(varargin{1},'ELLIPSE')
        u = varargin{1};
        u.a = varargin{2};
        u.b = varargin{3};
    end
elseif nargin==4
    u.dim = 2;
    u.cx = varargin{1};
    u.cy = varargin{2};
    u.cz = 0;
    u.a  = varargin{3};
    u.b  = varargin{4};
    u.nx = 0;
    u.ny = 0;
    u.nz = 1;
    u.vx = 1;
    u.vy = 0;
    u.indim = 2;
    
    u = class(u,'ELLIPSE',GEOMOBJECT(u.dim,u.indim));
elseif nargin==6
    u.dim = 2;
    u.cx = varargin{1};
    u.cy = varargin{2};
    u.cz = 0;
    u.a  = varargin{3};
    u.b  = varargin{4};
    u.nx = 0;
    u.ny = 0;
    u.nz = 1;
    u.vx = varargin{5};
    u.vy = varargin{6};
    u.indim = 2;
    
    u = class(u,'ELLIPSE',GEOMOBJECT(u.dim,u.indim));
elseif nargin==5
    u.dim = 2;
    u.cx = varargin{1};
    u.cy = varargin{2};
    u.cz = varargin{3};
    u.a  = varargin{4};
    u.b  = varargin{5};
    u.nx = 0;
    u.ny = 0;
    u.nz = 1;
    u.vx = 1;
    u.vy = 0;
    u.indim = 3;
    
    u = class(u,'ELLIPSE',GEOMOBJECT(u.dim,u.indim));
elseif nargin==7
    u.dim = 2;
    u.cx = varargin{1};
    u.cy = varargin{2};
    u.cz = varargin{3};
    u.a  = varargin{4};
    u.b  = varargin{5};
    u.nx = 0;
    u.ny = 0;
    u.nz = 1;
    u.vx = varargin{6};
    u.vy = varargin{7};
    u.indim = 3;
    
    u = class(u,'ELLIPSE',GEOMOBJECT(u.dim,u.indim));
elseif nargin==8
    u.dim = 2;
    u.cx = varargin{1};
    u.cy = varargin{2};
    u.cz = varargin{3};
    u.a  = varargin{4};
    u.b  = varargin{5};
    u.nx = varargin{6};
    u.ny = varargin{7};
    u.nz = varargin{8};
    u.vx = 1;
    u.vy = 0;
    u.indim = 3;
    
    u = class(u,'ELLIPSE',GEOMOBJECT(u.dim,u.indim));
elseif nargin==10
    u.dim = 2;
    u.cx = varargin{1};
    u.cy = varargin{2};
    u.cz = varargin{3};
    u.a  = varargin{4};
    u.b  = varargin{5};
    u.nx = varargin{6};
    u.ny = varargin{7};
    u.nz = varargin{8};
    u.vx = varargin{9};
    u.vy = varargin{10};
    u.indim = 3;
    
    u = class(u,'ELLIPSE',GEOMOBJECT(u.dim,u.indim));
% else
%     error('Wrong input arguments')
end
end
