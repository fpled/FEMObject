function u = CIRCLE(varargin)
% function C = CIRCLE(cx,cy,r)
% function C = CIRCLE(cx,cy,r,vx,vy)
% function C = CIRCLE(cx,cy,cz,r)
% function C = CIRCLE(cx,cy,cz,r,vx,vy)
% function C = CIRCLE(cx,cy,cz,r,nx,ny,nz)
% function C = CIRCLE(cx,cy,cz,r,nx,ny,nz,vx,vy)

if nargin==0
    u = CIRCLE(0,0,1);
elseif nargin==1
    if isa(varargin{1},'CIRCLE')
        u = varargin{1};
    end
elseif nargin==3
    u.dim = 2;
    u.cx = varargin{1};
    u.cy = varargin{2};
    u.cz = 0;
    u.r  = varargin{3};
    u.nx = 0;
    u.ny = 0;
    u.nz = 1;
    u.vx = 1;
    u.vy = 0;
    u.indim = 2;
    
    u = class(u,'CIRCLE',GEOMOBJECT(u.dim,u.indim));
elseif nargin==5
    u.dim = 2;
    u.cx = varargin{1};
    u.cy = varargin{2};
    u.cz = 0;
    u.r  = varargin{3};
    u.nx = 0;
    u.ny = 0;
    u.nz = 1;
    u.vx = varargin{4};
    u.vy = varargin{5};
    u.indim = 2;
    
    u = class(u,'CIRCLE',GEOMOBJECT(u.dim,u.indim));
elseif nargin==4
    u.dim = 2;
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
    
    u = class(u,'CIRCLE',GEOMOBJECT(u.dim,u.indim));
elseif nargin==6
    u.dim = 2;
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
    
    u = class(u,'CIRCLE',GEOMOBJECT(u.dim,u.indim));
elseif nargin==7
    u.dim = 2;
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
    
    u = class(u,'CIRCLE',GEOMOBJECT(u.dim,u.indim));
elseif nargin==9
    u.dim = 2;
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
    
    u = class(u,'CIRCLE',GEOMOBJECT(u.dim,u.indim));
% else
%     error('Wrong input arguments')
end
end
