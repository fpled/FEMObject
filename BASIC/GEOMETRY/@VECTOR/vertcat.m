function w = vertcat(u,v,varargin)
% function w = vertcat(u,v,varargin)

if isa(u,'VECTOR') && isa(v,'VECTOR')
    w = u;
    w.MYDOUBLEND = concat(u.MYDOUBLEND,v.MYDOUBLEND,3);
end

if length(varargin)>0
    w = vertcat(w,varargin{:});
end
