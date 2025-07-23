function lsx=eval(ls,x,varargin)

dim = (length(varargin)-2)/2;

repP = [1:dim];
repV = [dim+1:2*dim];
r=varargin{end-1};
r=MYDOUBLEND(reshape(full(r),[1,1,1,numel(r)]));
normtype = varargin{end};
for i=1:2*dim
    varargin{i}=MYDOUBLEND(reshape(full(varargin{i}),[1,1,1,numel(varargin{i})]));
end

P = POINT([varargin{repP}]);
V = VECTOR([varargin{repV}]');

lsx=distance(x,STRAIGHTLINE(P,V),normtype)-r;