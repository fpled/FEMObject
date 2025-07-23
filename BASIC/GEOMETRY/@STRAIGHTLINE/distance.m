function w = distance(D,v,varargin)
% function w = distance(D,v,varargin)

if isa(v,'POINT')
    
    P1 = getP(D,1);
    d = getV(D,1);
    P1M = v - P1;
    
    HM = P1M - d*dot(P1M,d);
    w = norm(HM,varargin{:})';
    
    if numel(w)==1
        w = double(w);
    end
    
end
