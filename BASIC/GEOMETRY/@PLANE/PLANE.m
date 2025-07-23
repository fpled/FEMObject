function u = PLANE(varargin)
%function D = PLANE(P1,P2,P3)
%       3 POINTS
%function D = PLANE(P1,V1)
%       1 POINT ET 1 VECTOR (normal au plan)
%function D = PLANE(P1,V1,V2)
%       1 POINT ET 2 VECTORS
%function D = PLANE(P1,P2,V1)
%       2 POINTS ET 1 VECTOR

if nargin==0
    P1 = POINT([0,0,0]);
    P2 = POINT([1,0,0]);
    P3 = POINT([1,1,0]);
    u = PLANE(P1,P2,P3);
elseif nargin==1
    if isa(varargin{1},'PLANE')
        u = varargin{1};
    end
else
    if isa(varargin{1},'double')
        varargin{1} = POINT(varargin{1});
    end
    if isa(varargin{2},'double')
        varargin{2} = POINT(varargin{2});
    end
    if nargin>=3 && isa(varargin{3},'double')
        varargin{3} = POINT(varargin{3});
    end
    
    [V,posV,nV] = getclassin('VECTOR',varargin);
    [P,posP,nP] = getclassin('POINT',varargin);
    
    if ~isa(V,'cell')
        V = {V};
    end
    if ~isa(P,'cell')
        P = {P};
    end
    
    if (nV==2 && nP==1)
        u.P{1} = P{1};
        u.V{2} = normalize(V{1});
        u.V{3} = normalize(V{2});
        u.P{2} = u.P{1}+u.V{2};
        u.P{3} = u.P{1}+u.V{3};
        
    elseif (nV==1 && nP==1)
        v1 = V{1};
        u.P{1} = P{1};
        [v2,v3] = planortho(v1);
        u.V{2} = v2;
        u.V{3} = v3;
        u.P{2} = u.P{1}+u.V{2};
        u.P{3} = u.P{1}+u.V{3};
        
    elseif (nV==1 && nP==2)
        u.P{1} = P{1};
        u.P{2} = P{2};
        u.V{2} = normalize(u.P{2} - u.P{1}) ;
        u.V{3} = normalize(V{1});
        u.P{3} = u.P{1}+u.V{3}*distance(u.P{1},u.P{2});
        
    elseif (nV==0 && nP==3)
        u.P{1} = P{1};
        u.P{2} = P{2};
        u.P{3} = P{3};
        u.V{2} = normalize(u.P{2} - u.P{1}) ;
        u.V{3} = normalize(u.P{3} - u.P{1});
        
    else
        help PLANE
        error('Wrong argument')
    end
    
    if norm(cross(u.V{2},u.V{3}))<10*eps
        error('Collinear vectors')
    else
        u.V{1} = cross(u.V{2},u.V{3});
    end
    
    u.dim = 2;
    u.indim = 3;
    
    u = class(u,'PLANE',GEOMOBJECT(u.dim,u.indim));
end
