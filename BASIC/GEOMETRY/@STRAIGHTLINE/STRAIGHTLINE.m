function u = STRAIGHTLINE(varargin)
% function D = STRAIGHTLINE(P1,P2)
%       2 POINTS
% function D = STRAIGHTLINE(P,V)
%       1 POINT ET 1 VECTOR

if nargin==0
    u.P = cell(1,2);
    u.dim = [];
    u.indim = [];
    u.V = cell(1,0);
    u = class(u,'STRAIGHTLINE',GEOMOBJECT());
elseif nargin==1
    if isa(varargin{1},'STRAIGHTLINE')
        u = varargin{1};
    end
elseif nargin==2
    if isa(varargin{1},'double')
        varargin{1} = POINT(varargin{1});
    end
    if isa(varargin{2},'double')
        varargin{2} = POINT(varargin{2});
    end
    if isa(varargin{1},'POINT') && isa(varargin{2},'POINT')
        P1 = varargin{1};
        P2 = varargin{2};
        L = distance(P1,P2);
        v = (P2-P1)/L;
        
    elseif isa(varargin{1},'POINT') && isa(varargin{2},'VECTOR')
        P1 = varargin{1};
        v = normalize(varargin{2});
        P2 = P1 + v;
        
    else
        help STRAIGHTLINE
        error('rentrer les bons arguments')
    end
    u.P{1} = P1;
    u.P{2} = P2;
    u.V{1} = v;
    u.dim = 1;
    u.indim = getindim(u.P{1});
    
    u = class(u,'STRAIGHTLINE',GEOMOBJECT(u.dim,u.indim));
        
end
