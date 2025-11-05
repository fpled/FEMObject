function [se,B] = sigma(mat,elem,xnode,xgauss,qe,varargin)
% function [se,B] = sigma(mat,elem,xnode,xgauss,qe,varargin)

B = calc_B(elem,xnode,xgauss);

if isparam(mat,'d')
    model = getparam(mat,'PFM'); % phase field model
    
    switch lower(model)
        case 'bourdin'
            Dp = calc_opmat(mat,elem,xnode,xgauss);
            Dm = 0;
        case 'amor'
            se = B*qe; % strain field
            [Dp,Dm] = calc_opmat_Amor(mat,elem,xnode,xgauss,se);
        case 'miehe'
            se = B*qe; % strain field
            [Dp,Dm] = calc_opmat_Miehe(mat,elem,xnode,xgauss,se);
        case {'heamor','hefreddi'}
            se = B*qe; % strain field
            [Dp,Dm] = calc_opmat_He(mat,elem,xnode,xgauss,se);
        case 'zhang'
            se = B*qe; % strain field
            [Dp,Dm] = calc_opmat_Zhang(mat,elem,xnode,xgauss,se);
        case 'spectral'
            se = B*qe; % strain field
            [Dp,Dm] = calc_opmat_spectral(mat,elem,xnode,xgauss,se);
        case 'doublespectral'
            [Dp,Dm] = calc_opmat_doublespectral(mat,elem,xnode,xgauss);
        otherwise
            error(['Wrong phase field model ' model])
    end
    
    if ischarin('positive',varargin)
        D = Dp;
    elseif ischarin('negative',varargin)
        D = Dm;
    else
        d = evalparam(mat,'d',elem,xnode,xgauss); % phase field
        if isparam(mat,'h')
            h = evalparam(mat,'h',elem,xnode,xgauss); % healing field
        end
        g = getparam(mat,'g'); % energetic degradation function
        k = evalparam(mat,'k',elem,xnode,xgauss); % small artificial residual stiffness
        
        N = calc_N(elem,xnode,xgauss,'nbddlpernode',1);
        de = localize(elem,d,'scalar');
        de = N*de;
        if isparam(mat,'h')
            he = localize(elem,h,'scalar');
            he = N*he;
            D = (g(de,he)+k)*Dp + Dm;
        else
            D = (g(de)+k)*Dp + Dm;
        end
    end
else
    D = calc_opmat(mat,elem,xnode,xgauss);
end

se = D*(B*qe);

if ischarin('local',varargin)
    switch getdim(elem)
        case 1
            S = evalparam(mat,'S',elem,xnode,xgauss); % cross-section area
            se = se/S;
        case 2
            e = evalparam(mat,'DIM3',elem,xnode,xgauss); % thickness
            se = se/e;
    end
end
