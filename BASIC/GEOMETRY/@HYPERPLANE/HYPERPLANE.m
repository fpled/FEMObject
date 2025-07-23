function u = HYPERPLANE(varargin)
%function D = HYPERPLANE(P1,V1)
%       1 POINT ET 1 VECTOR (normal au plan)

switch nargin
    case 1
        if isa(varargin{1},'HYPERPLANE')
            u = varargin{1};
        end
    otherwise
        
        [V,nV] = getclassinvarargin('VECTOR',varargin);
        [P,nP] =  getclassinvarargin('POINT',varargin);
        if ~isa(V,'cell')
            V = {V};
        end
        if ~isa(P,'cell')
            P = {P};
        end
        if nV==1 && nP==1
            switch getindim(P{1})
                case 1
                    u.P{1} = P{1};
                    u.V{1} = V{1};
                    u.dim = 0;   
                case 2
                    u.P{1} = P{1};
                    u.V{1} = V{1};
                    u.V{2} = rot2D(V{1},pi/2);
                    u.dim = 1;
                case 3
                    u.P{1} = P{1};
                    u.V{1} = V{1};
                    [u.V{2},u.V{3}] = planortho(V{1});
                    u.dim = 2;
            end
        else
            error('Wrong input arguments')
        end
        
        u.indim = getindim(u.P{1});
        u = class(u,'HYPERPLANE',GEOMOBJECT(u.dim,u.indim));
        
end
