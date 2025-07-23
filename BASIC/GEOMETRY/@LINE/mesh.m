function S = mesh(B,varargin)
% function S = mesh(B,n)
% function S = mesh(B,n,'indim',indim...)
% n : numbre of elements
% indim : space dimension (B.indim by default)
% optional arguments for addelem : material, param, option

indim = getcharin('indim',varargin,B.indim);
varargin = delcharin('indim',varargin);

switch indim
    case 1
        S = MODEL('UNID');
    case 2
        S = MODEL('PLANE');
    case 3
        S = MODEL('TRID');
    otherwise
        error('Wrong space dimension')
end

n1 = varargin{1};

v = B.P{2}-B.P{1};
L = distance(B.P{1},B.P{2}) ;
v = normalize(v);
X = linspace(0,L,n1+1);

X = reshape(X,[1 1 numel(X)]);

node = NODE(B.P{1}+v*X,1:n1+1);
elem = [1:n1;2:n1+1]';

S = addnode(S,node);
S = addelem(S,'SEG2',elem,varargin{:});
