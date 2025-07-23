function S = mesh(B,varargin)
% function S = mesh(B,n1,n2)
% function S = mesh(B,n1,n2,'indim',indim,...)
% ni : numbre of elements in dimension i
% indim : space dimension (B.indim by default)
% optional arguments for addelem : material, param, option

indim = getcharin('indim',varargin,B.indim);
varargin = delcharin('indim',varargin);

switch indim
    case 2
        S = MODEL('PLANE');
    case 3
        S = MODEL('TRID');
    otherwise
        error('Wrong space dimension')
end

n1 = varargin{1};
n2 = varargin{2};

P = vertcat(B.P{:});
v1 = P(2)-P(1);
v2 = P(4)-P(1);
node = NODE(P,1:4);
elem = QUA4(node,1,1:4);

[xi,eta] = meshgrid(linspace(-1,1,n1+1),linspace(-1,1,n2+1));

xnode = calc_x(elem,getcoord(node),[xi(:),eta(:)]);

node = NODE(xnode,1:length(xi(:)));

elem = zeros(n1*n2,4);
n = 0;
for j=1:n1
    i = (1:n2)';
    connec = [(j-1)*(n2+1)+i,...
        (j)*(n2+1)+i,...
        (j)*(n2+1)+i+1,...
        (j-1)*(n2+1)+i+1]    ;
    elem(n+(1:n2),:) = connec;
    n = n+n2;
end

S = addnode(S,node);
S = addelem(S,'QUA4',elem,varargin{:});
