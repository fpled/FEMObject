function X2 = transfer(D1,D2,X1,dim)
% function X2 = transfer(D1,D2,X1)
% transformation between D1 and D2: map the points X1 into points X2
% D1 and D2 : DOMAIN of dimension d
% X1 a n-by-d array containing coordinates of points in the domain D1
% X2 a n-by-d array containing the corresponding coordinates in D2
%
% function X2 = transfer(D1,D2,X1,dim)
% only the coordinates dim are given for points

a1 = double(getcoord(D1.P1));
b1 = double(getcoord(D1.P2));
a2 = double(getcoord(D2.P1));
b2 = double(getcoord(D2.P2));

if nargin<4 || isempty(dim)
    dim = 1:getdim(D1);
end

X2 = zeros(size(X1));
for k=1:length(dim)
    i = dim(k);
    Y = (X1(:,k)-a1(i))/(b1(i)-a1(i));
    X2(:,k) = a2(i) + (b2(i)-a2(i))*Y;
end
