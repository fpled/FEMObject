function w = distance(u,v,varargin)
% function w = distance(u,v,normtype)
% distance de u a v
% normtype : 1, 2 ou Inf (type de norme) 2 par defaut
% u et v peuvent etre des POINT ou des GEOMOBJECT
%
% See also POINT/distance, POINT/minus, POINT/mtimes, POINT/mrdivide, POINT/ne, POINT/eq,
% POINT/plus, POINT/uminus, POINT/norm, VECTOR/ne,
% PLANE/distance, HYPERPLANE/distance, STRAIGHTLINE/distance

if isa(v,'POINT')
    w = norm(POINT(u)-v,varargin{:})';
    if numel(w)==1
        w = double(w);
    end
    
elseif isa(v,'GEOMOBJECT')
    w = distance(v,u,varargin{:});
    
else
    error('pas programme')
end

