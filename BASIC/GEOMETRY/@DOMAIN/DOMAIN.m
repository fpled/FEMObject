function u = DOMAIN(varargin)
% function D = DOMAIN(dim,P1,P2)
% Pi : arguments de type double (2 points extremes)
% dim : dimension du domaine (peut etre un domaine 1D, 2D ou 3D)

if nargin==0
    u = DOMAIN(1);
elseif nargin==1
    if isa(varargin{1},'DOMAIN')
        u = varargin{1};
    elseif isa(varargin{1},'double')
        P1 = zeros(1,varargin{1});
        P2 = ones(1,varargin{1});
        u = DOMAIN(varargin{1},P1,P2) ;
    end
elseif nargin==3
    u.dim = varargin{1};
    u.P1 = POINT(varargin{2});
    u.P2 = POINT(varargin{3});
    u.indim = getindim(u.P1);
    u.diameter = distance(u.P1,u.P2);
    if ~(u.indim==u.dim)
        if u.dim==1
            error('Use LINE')
        elseif u.dim==2
            error('Use QUADRANGLE')
        end
    end
    u = class(u,'DOMAIN',GEOMOBJECT(u.dim,u.indim));
end
