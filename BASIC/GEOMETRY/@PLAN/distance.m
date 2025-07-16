function w = distance(u,P,varargin)
% function w = distance(u,P,varargin)

if isa(P,'POINT')
    
    P1P = P-u.P{1};
    w = abs(dot(P1P,u.V{1}));
    
    if isscalar(w)
        w = double(w);
    end
    
end
