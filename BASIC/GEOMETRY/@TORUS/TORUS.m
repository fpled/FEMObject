function u = TORUS(varargin)
% function T = TORUS(T)
% function T = TORUS(T,r2)
% function T = TORUS(C,r2)
% function T = TORUS(T,r2,angle)
% function T = TORUS(C,r2,angle)
% function T = TORUS(cx,cy,cz,r1,r2)
% function T = TORUS(cx,cy,cz,r1,r2,angle)
% function T = TORUS(cx,cy,cz,r1,r2,vx,vy)
% function T = TORUS(cx,cy,cz,r1,r2,nx,ny,nz)
% function T = TORUS(cx,cy,cz,r1,r2,nx,ny,nz,angle)
% function T = TORUS(cx,cy,cz,r1,r2,nx,ny,nz,vx,vy)
% function T = TORUS(cx,cy,cz,r1,r2,nx,ny,nz,vx,vy,angle)

if nargin==0
    u = TORUS(0,0,1,1);
elseif nargin==1
    if isa(varargin{1},'TORUS')
        u = varargin{1};
    end
elseif nargin==2
    if isa(varargin{1},'TORUS')
        u = varargin{1};
        u.r2 = varargin{2};
    elseif isa(varargin{1},'CIRCLE')
        C = varargin{1};
        c = getc(C);
        r = getr(C);
        n = getn(C);
        v = getv(C);
        u.dim = 3;
        u.cx = c(1);
        u.cy = c(2);
        u.cz = c(3);
        u.r1 = r;
        u.r2 = varargin{2};
        u.nx = n(1);
        u.ny = n(2);
        u.nz = n(3);
        u.vx = v(1);
        u.vy = v(2);
        u.angle = 2*pi;
        u.indim = 3;
        
        u = class(u,'TORUS',GEOMOBJECT(u.dim,u.indim));
    end
elseif nargin==3
    if isa(varargin{1},'TORUS')
        u = varargin{1};
        u.r2 = varargin{2};
        u.angle = varargin{2};
    elseif isa(varargin{1},'CIRCLE')
        C = varargin{1};
        c = getc(C);
        r = getr(C);
        n = getn(C);
        v = getv(C);
        u.dim = 3;
        u.cx = c(1);
        u.cy = c(2);
        u.cz = c(3);
        u.r1 = r;
        u.r2 = varargin{2};
        u.nx = n(1);
        u.ny = n(2);
        u.nz = n(3);
        u.vx = v(1);
        u.vy = v(2);
        u.angle = varargin{3};
        u.indim = 3;
        
        u = class(u,'TORUS',GEOMOBJECT(u.dim,u.indim));
    end
elseif nargin==5
    u.dim = 3;
    u.cx = varargin{1};
    u.cy = varargin{2};
    u.cz = varargin{3};
    u.r1 = varargin{4};
    u.r2 = varargin{5};
    u.nx = 0;
    u.ny = 0;
    u.nz = 1;
    u.vx = 1;
    u.vy = 0;
    u.angle = 2*pi;
    u.indim = 3;
    
    u = class(u,'TORUS',GEOMOBJECT(u.dim,u.indim));
elseif nargin==6
    u.dim = 3;
    u.cx = varargin{1};
    u.cy = varargin{2};
    u.cz = varargin{3};
    u.r1 = varargin{4};
    u.r2 = varargin{5};
    u.nx = 0;
    u.ny = 0;
    u.nz = 1;
    u.vx = 1;
    u.vy = 0;
    u.angle = varargin{6};
    u.indim = 3;
    
    u = class(u,'TORUS',GEOMOBJECT(u.dim,u.indim));
elseif nargin==7
    u.dim = 3;
    u.cx = varargin{1};
    u.cy = varargin{2};
    u.cz = varargin{3};
    u.r1 = varargin{4};
    u.r2 = varargin{5};
    u.nx = 0;
    u.ny = 0;
    u.nz = 1;
    u.vx = varargin{6};
    u.vy = varargin{7};
    u.angle = 2*pi;
    u.indim = 3;

    u = class(u,'TORUS',GEOMOBJECT(u.dim,u.indim));
elseif nargin==8
    u.dim = 3;
    u.cx = varargin{1};
    u.cy = varargin{2};
    u.cz = varargin{3};
    u.r1 = varargin{4};
    u.r2 = varargin{5};
    u.nx = varargin{6};
    u.ny = varargin{7};
    u.nz = varargin{8};
    u.vx = 1;
    u.vy = 0;
    u.angle = 2*pi;
    u.indim = 3;
    
    u = class(u,'TORUS',GEOMOBJECT(u.dim,u.indim));
elseif nargin==9
    u.dim = 3;
    u.cx = varargin{1};
    u.cy = varargin{2};
    u.cz = varargin{3};
    u.r1 = varargin{4};
    u.r2 = varargin{5};
    u.nx = varargin{6};
    u.ny = varargin{7};
    u.nz = varargin{8};
    u.vx = 1;
    u.vy = 0;
    u.angle = varargin{9};
    u.indim = 3;
    
    u = class(u,'TORUS',GEOMOBJECT(u.dim,u.indim));
elseif nargin==10
    u.dim = 3;
    u.cx = varargin{1};
    u.cy = varargin{2};
    u.cz = varargin{3};
    u.r1 = varargin{4};
    u.r2 = varargin{5};
    u.nx = varargin{6};
    u.ny = varargin{7};
    u.nz = varargin{8};
    u.vx = varargin{9};
    u.vy = varargin{10};
    u.angle = 2*pi;
    u.indim = 3;
    
    u = class(u,'TORUS',GEOMOBJECT(u.dim,u.indim));
elseif nargin==11
    u.dim = 3;
    u.cx = varargin{1};
    u.cy = varargin{2};
    u.cz = varargin{3};
    u.r1 = varargin{4};
    u.r2 = varargin{5};
    u.nx = varargin{6};
    u.ny = varargin{7};
    u.nz = varargin{8};
    u.vx = varargin{9};
    u.vy = varargin{10};
    u.angle = varargin{11};
    u.indim = 3;
    
    u = class(u,'TORUS',GEOMOBJECT(u.dim,u.indim));
% else
%     error('Wrong input arguments')
end
end
