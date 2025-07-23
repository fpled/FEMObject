function lsx=eval(ls,x,varargin)

dim = (length(varargin))/2;

repP = 1:dim;
repV = dim+[1:dim];
for i=1:2*dim
    varargin{i}=MYDOUBLEND(reshape(full(varargin{i}),[1,1,1,numel(varargin{i})]));
end

P = POINT([varargin{repP}]);
V = VECTOR([varargin{repV}]');

lsx=dot(V,POINT(x) - P);
