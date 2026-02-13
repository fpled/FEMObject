function D = calc_opmasse(mat,elem,xnode,xgauss)
% function D = calc_opmasse(mat,elem,xnode,xgauss)

rho = evalparam(mat,'RHO',elem,xnode,xgauss); % mass density

switch getdim(elem)
    case 1
        S = evalparam(mat,'S',elem,xnode,xgauss); % cross-section area
        D = rho*S;
    case 2
        e = evalparam(mat,'DIM3',elem,xnode,xgauss); % thickness
        D = rho*e;
    case 3
        D = rho;
end
