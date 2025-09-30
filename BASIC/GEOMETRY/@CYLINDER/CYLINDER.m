function u = CYLINDER(varargin)
% function C = CYLINDER(C)
% function C = CYLINDER(C,h)
% function C = CYLINDER(C,h,angle)
% function C = CYLINDER(cx,cy,cz,r,h)
% function C = CYLINDER(cx,cy,cz,r,h,angle)
% function C = CYLINDER(cx,cy,cz,r,h,vx,vy)
% function C = CYLINDER(cx,cy,cz,r,h,nx,ny,nz)
% function C = CYLINDER(cx,cy,cz,r,h,nx,ny,nz,angle)
% function C = CYLINDER(cx,cy,cz,r,h,nx,ny,nz,vx,vy)
% function C = CYLINDER(cx,cy,cz,r,h,nx,ny,nz,vx,vy,angle)

if nargin==0
    u = CYLINDER(0,0,1,1);
elseif nargin==1
    if isa(varargin{1},'CYLINDER')
        u = varargin{1};
    end
elseif nargin==2
    if isa(varargin{1},'CYLINDER')
        u = varargin{1};
        u.h = varargin{2};
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
        u.r  = r;
        u.h  = varargin{2};
        u.nx = n(1);
        u.ny = n(2);
        u.nz = n(3);
        u.vx = v(1);
        u.vy = v(2);
        u.angle = 2*pi;
        u.indim = 3;
        
        u = class(u,'CYLINDER',GEOMOBJECT(u.dim,u.indim));
    end
elseif nargin==3
    if isa(varargin{1},'CYLINDER')
        u = varargin{1};
        u.r = varargin{2};
        u.h = varargin{3};
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
        u.r  = r;
        u.h  = varargin{2};
        u.nx = n(1);
        u.ny = n(2);
        u.nz = n(3);
        u.vx = v(1);
        u.vy = v(2);
        u.angle = varargin{3};
        u.indim = 3;
        
        u = class(u,'CYLINDER',GEOMOBJECT(u.dim,u.indim));
    end
elseif nargin==5
    u.dim = 3;
    u.cx = varargin{1};
    u.cy = varargin{2};
    u.cz = varargin{3};
    u.r  = varargin{4};
    u.h  = varargin{5};
    u.nx = 0;
    u.ny = 0;
    u.nz = 1;
    u.vx = 1;
    u.vy = 0;
    u.angle = 2*pi;
    u.indim = 3;
    
    u = class(u,'CYLINDER',GEOMOBJECT(u.dim,u.indim));
elseif nargin==6
    u.dim = 3;
    u.cx = varargin{1};
    u.cy = varargin{2};
    u.cz = varargin{3};
    u.r  = varargin{4};
    u.h  = varargin{5};
    u.nx = 0;
    u.ny = 0;
    u.nz = 1;
    u.vx = 1;
    u.vy = 0;
    u.angle = varargin{6};
    u.indim = 3;
    
    u = class(u,'CYLINDER',GEOMOBJECT(u.dim,u.indim));
elseif nargin==7
    u.dim = 3;
    u.cx = varargin{1};
    u.cy = varargin{2};
    u.cz = varargin{3};
    u.r  = varargin{4};
    u.h  = varargin{5};
    u.nx = 0;
    u.ny = 0;
    u.nz = 1;
    u.vx = varargin{6};
    u.vy = varargin{7};
    u.angle = 2*pi;
    u.indim = 3;
    
    u = class(u,'CYLINDER',GEOMOBJECT(u.dim,u.indim));
elseif nargin==8
    u.dim = 3;
    u.cx = varargin{1};
    u.cy = varargin{2};
    u.cz = varargin{3};
    u.r  = varargin{4};
    u.h  = varargin{5};
    u.nx = varargin{6};
    u.ny = varargin{7};
    u.nz = varargin{8};
    u.vx = 1;
    u.vy = 0;
    u.angle = 2*pi;
    u.indim = 3;
    
    u = class(u,'CYLINDER',GEOMOBJECT(u.dim,u.indim));
elseif nargin==9
    u.dim = 3;
    u.cx = varargin{1};
    u.cy = varargin{2};
    u.cz = varargin{3};
    u.r  = varargin{4};
    u.h  = varargin{5};
    u.nx = varargin{6};
    u.ny = varargin{7};
    u.nz = varargin{8};
    u.vx = 1;
    u.vy = 0;
    u.angle = varargin{9};
    u.indim = 3;
    
    u = class(u,'CYLINDER',GEOMOBJECT(u.dim,u.indim));
elseif nargin==10
    u.dim = 3;
    u.cx = varargin{1};
    u.cy = varargin{2};
    u.cz = varargin{3};
    u.r  = varargin{4};
    u.h  = varargin{5};
    u.nx = varargin{6};
    u.ny = varargin{7};
    u.nz = varargin{8};
    u.vx = varargin{9};
    u.vy = varargin{10};
    u.angle = 2*pi;
    u.indim = 3;
    
    u = class(u,'CYLINDER',GEOMOBJECT(u.dim,u.indim));
elseif nargin==11
    u.dim = 3;
    u.cx = varargin{1};
    u.cy = varargin{2};
    u.cz = varargin{3};
    u.r  = varargin{4};
    u.h  = varargin{5};
    u.nx = varargin{6};
    u.ny = varargin{7};
    u.nz = varargin{8};
    u.vx = varargin{9};
    u.vy = varargin{10};
    u.angle = varargin{11};
    u.indim = 3;
    
    u = class(u,'CYLINDER',GEOMOBJECT(u.dim,u.indim));
% else
%     error('Wrong input arguments')
end
end
