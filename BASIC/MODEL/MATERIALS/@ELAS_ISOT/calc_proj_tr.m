function [Rp,Rm] = calc_proj_tr(mat,elem,xnode,xgauss,se)
% function [Rp,Rm] = calc_proj_tr(mat,elem,xnode,xgauss,se)

switch getdim(elem)
    case 1
        trSe = se(1);
    case 2
        trSe = se(1)+se(2);
    case 3
        trSe = se(1)+se(2)+se(3);
end

Rp = heaviside(trSe);  % Rp = (1+sign(trSe))/2;
Rm = heaviside(-trSe); % Rm = (1+sign(-trSe))/2;

end
