function u = LINE(varargin)
%function L = LINE(P1,P2)
% Pi : arguments de type POINT

if nargin==0
    u.P = cell(1,2);
    u.dim = [];
    u.indim = [];
    u.L = [];
    u.V = cell(1,0);
    u = class(u,'LINE',GEOMOBJECT());
elseif nargin==1
    if isa(varargin{1},'LINE')
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
        u.P{1} = varargin{1};
        u.P{2} = varargin{2};
        u.dim = 1;
        u.indim = getindim(u.P{1});

        u.L = distance(u.P{1},u.P{2});
        u.V{1} = (u.P{2}-u.P{1})/u.L;
        u = class(u,'LINE',GEOMOBJECT(u.dim,u.indim));
    end

end
