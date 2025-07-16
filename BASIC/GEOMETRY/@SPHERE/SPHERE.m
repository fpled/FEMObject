function u = SPHERE(varargin)
% function S = SPHERE(C)
% function S = SPHERE(C,r)
% function S = SPHERE(cx,cy,cz,r)
% function S = SPHERE(cx,cy,cz,r,vx,vy)
% function S = SPHERE(cx,cy,cz,r,nx,ny,nz)
% function S = SPHERE(cx,cy,cz,r,nx,ny,nz,vx,vy)

if nargin==0
    u = SPHERE(0,0,0,1);
elseif nargin==1
    if isa(varargin{1},'SPHERE')
        u = varargin{1};
    elseif isa(varargin{1},'CIRCLE')
        C = varargin{1};
        c = getc(C);
        n = getn(C);
        v = getv(C);
        r = getr(C);
        u.dim = 3;
        u.cx = c(1);
        u.cy = c(2);
        u.cz = c(3);
        u.r  = r;
        u.nx = n(1);
        u.ny = n(2);
        u.nz = n(3);
        u.vx = v(1);
        u.vy = v(2);
        u.indim = 3;
        
        u = class(u,'SPHERE',GEOMOBJECT(u.dim,u.indim));
    end
elseif nargin==2
    if isa(varargin{1},'CIRCLE')
        C = varargin{1};
        c = getc(C);
        n = getn(C);
        v = getv(C);
        % r = getr(C);
        u.dim = 3;
        u.cx = c(1);
        u.cy = c(2);
        u.cz = c(3);
        u.r  = varargin{2};
        u.nx = n(1);
        u.ny = n(2);
        u.nz = n(3);
        u.vx = v(1);
        u.vy = v(2);
        u.indim = 3;
        
        u = class(u,'SPHERE',GEOMOBJECT(u.dim,u.indim));
    end
elseif nargin==4
    u.dim = 3;
    u.cx = varargin{1};
    u.cy = varargin{2};
    u.cz = varargin{3};
    u.r  = varargin{4};
    u.nx = 0;
    u.ny = 0;
    u.nz = 1;
    u.vx = 1;
    u.vy = 0;
    u.indim = 3;
    
    u = class(u,'SPHERE',GEOMOBJECT(u.dim,u.indim));
elseif nargin==6
    u.dim = 3;
    u.cx = varargin{1};
    u.cy = varargin{2};
    u.cz = varargin{3};
    u.r  = varargin{4};
    u.nx = 0;
    u.ny = 0;
    u.nz = 1;
    u.vx = varargin{5};
    u.vy = varargin{6};
    u.indim = 3;
    
    u = class(u,'SPHERE',GEOMOBJECT(u.dim,u.indim));
elseif nargin==7
    u.dim = 3;
    u.cx = varargin{1};
    u.cy = varargin{2};
    u.cz = varargin{3};
    u.r  = varargin{4};
    u.nx = varargin{5};
    u.ny = varargin{6};
    u.nz = varargin{7};
    u.vx = 1;
    u.vy = 0;
    u.indim = 3;
    
    u = class(u,'SPHERE',GEOMOBJECT(u.dim,u.indim));
elseif nargin==9
    u.dim = 3;
    u.cx = varargin{1};
    u.cy = varargin{2};
    u.cz = varargin{3};
    u.r  = varargin{4};
    u.nx = varargin{5};
    u.ny = varargin{6};
    u.nz = varargin{7};
    u.vx = varargin{8};
    u.vy = varargin{9};
    u.indim = 3;
    
    u = class(u,'SPHERE',GEOMOBJECT(u.dim,u.indim));
% else
%     error('Wrong input arguments')
end
end
