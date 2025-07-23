function z = end(u,k,n)
% function z = end(u,k,n)

switch n
    case 1
        z = prod(sizeND(u));
    case 2
        z = size(u,k+2);
end
