function DN = getDN(elem,xi)
% function DN = getDN(elem,xi)

if nargin==2
    s = sizeND(xi);
    xi = double(xi);
    DN = zeros([3,8,sizeND(xi)]);
    
    fact1 = repmat([-1,1,1,-1,-1,1,1,-1],[1,1,prod(s)]);
    fact2 = repmat([-1,-1,1,1,-1,-1,1,1],[1,1,prod(s)]);
    fact3 = repmat([-1,-1,-1,-1,1,1,1,1],[1,1,prod(s)]);
    xi1 = repmat(xi(1,1,:),[1,8,1]);
    xi2 = repmat(xi(1,2,:),[1,8,1]);
    xi3 = repmat(xi(1,3,:),[1,8,1]);
    
    DN(1,:,:) = 1/8*fact1.*(1+xi2.*fact2).*(1+xi3.*fact3);
    DN(2,:,:) = 1/8*fact2.*(1+xi1.*fact1).*(1+xi3.*fact3);
    DN(3,:,:) = 1/8*fact3.*(1+xi1.*fact1).*(1+xi2.*fact2);
else
    % DN = inline('1/8*[[-1,1,1,-1,-1,1,1,-1].*(1+xi(2)*[-1,-1,1,1,-1,-1,1,1]).*(1+xi(3)*[-1,-1,-1,-1,1,1,1,1]);(1+xi(1)*[-1,1,1,-1,-1,1,1,-1]).*[-1,-1,1,1,-1,-1,1,1].*(1+xi(3)*[-1,-1,-1,-1,1,1,1,1]);(1+xi(1)*[-1,1,1,-1,-1,1,1,-1]).*(1+xi(2)*[-1,-1,1,1,-1,-1,1,1]).*[-1,-1,-1,-1,1,1,1,1]]','xi');
    DN = @(xi) 1/8*[[-1,1,1,-1,-1,1,1,-1].*(1+xi(2)*[-1,-1,1,1,-1,-1,1,1]).*(1+xi(3)*[-1,-1,-1,-1,1,1,1,1]);(1+xi(1)*[-1,1,1,-1,-1,1,1,-1]).*[-1,-1,1,1,-1,-1,1,1].*(1+xi(3)*[-1,-1,-1,-1,1,1,1,1]);(1+xi(1)*[-1,1,1,-1,-1,1,1,-1]).*(1+xi(2)*[-1,-1,1,1,-1,-1,1,1]).*[-1,-1,-1,-1,1,1,1,1]];
end