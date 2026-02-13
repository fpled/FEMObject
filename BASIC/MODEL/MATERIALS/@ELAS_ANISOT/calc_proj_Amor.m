function [Pp,Pm] = calc_proj_Amor(mat,elem,xnode,xgauss,se,varargin)
% function [Pp,Pm] = calc_proj_Amor(mat,elem,xnode,xgauss,se)

dim = getdim(elem);
split = getparam(mat,'PFS'); % phase field split

P = calc_proj_notation(elem);
if strcmpi(split,'stress')
    P2 = P*P;
    se = P2*se; % stress tensor converted as strain tensor in Voigt notation
end

%% Computation of projectors of strain/stress tensor onto volumetric/spherical/hydrostatic expansion (or tensile bulk) and deviatoric parts
switch dim
    case 1
        v = 1;
        Id = 1;
    case 2
        if isaxi(elem)
            v = [1,1,1,0]';
            Id = diag([1,1,1,1/2]);
        else
            v = [1,1,0]';
            Id = diag([1,1,1/2]);
        end
    case 3
        v = [1,1,1,0,0,0]';
        Id = diag([1,1,1,1/2,1/2,1/2]);
end
I = v*v';
% n = length(v);
% Id = P'\eye(n)/P;

[~,Rm] = calc_proj_tr(mat,elem,xnode,xgauss,se);
Pm = (1/dim)*Rm*I;
Pp = Id - Pm;

%% Check orthogonality condition and strain/stress decomposition
if ischarin('check',varargin)
    tol = 1e-12;
    sep = P*(Pp*se); % positive part of strain/stress tensor in Kelvin-Mandel notation
    sem = P*(Pm*se); % negative part of strain/stress tensor in Kelvin-Mandel notation
    se = P\se; % strain/stress tensor in Kelvin-Mandel notation
    orthpm = max(abs(sep'*sem)/abs(se'*se),[],'all'); if orthpm>tol, orthpm, end
    orthmp = max(abs(sem'*sep)/abs(se'*se),[],'all'); if orthmp>tol, orthmp, end
    decomp = max(norm(se - (sep+sem))/norm(se),[],'all'); if decomp>tol, decomp, end
end

end
