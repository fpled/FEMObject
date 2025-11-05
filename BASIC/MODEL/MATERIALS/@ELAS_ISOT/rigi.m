function ke = rigi(mat,elem,xnode,xgauss,varargin)
% function ke = rigi(mat,elem,xnode,xgauss,varargin)

B = calc_B(elem,xnode,xgauss);

if isparam(mat,'d')
    model = getparam(mat,'PFM'); % phase field model
    
    switch lower(model)
        case 'bourdin'
            Dp = calc_opmat(mat,elem,xnode,xgauss);
            Dm = 0;
        case 'amor'
            u = evalparam(mat,'u',elem,xnode,xgauss); % displacement field
            ue = localize(elem,u);
            se = B*ue; % strain field
            [Dp,Dm] = calc_opmat_Amor(mat,elem,xnode,xgauss,se);
        case 'miehe'
            u = evalparam(mat,'u',elem,xnode,xgauss); % displacement field
            ue = localize(elem,u);
            se = B*ue; % strain field
            [Dp,Dm] = calc_opmat_Miehe(mat,elem,xnode,xgauss,se);
        case {'heamor','hefreddi'}
            u = evalparam(mat,'u',elem,xnode,xgauss); % displacement field
            ue = localize(elem,u);
            se = B*ue; % strain field
            [Dp,Dm] = calc_opmat_He(mat,elem,xnode,xgauss,se);
        case 'zhang'
            u = evalparam(mat,'u',elem,xnode,xgauss); % displacement field
            ue = localize(elem,u);
            se = B*ue; % strain field
            [Dp,Dm] = calc_opmat_Zhang(mat,elem,xnode,xgauss,se);
        case 'spectral'
            u = evalparam(mat,'u',elem,xnode,xgauss); % displacement field
            ue = localize(elem,u);
            se = B*ue; % strain field
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

ke = B'*D*B;
