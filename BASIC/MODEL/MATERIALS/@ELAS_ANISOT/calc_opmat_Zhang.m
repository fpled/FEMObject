function [Dp,Dm] = calc_opmat_Zhang(mat,elem,xnode,xgauss,se,varargin)
% function [Dp,Dm] = calc_opmat_Zhang(mat,elem,xnode,xgauss,se)

if nargin<=2
    xnode = [];
    xgauss = [];
    se = [];
end

dim = getdim(elem);
split = getparam(mat,'PFS'); % phase field split

D = calc_opmat(mat,elem,xnode,xgauss); % stiffness operator in Voigt notation

if dim==1
    if strcmpi(split,'stress')
        C = inv(D); % compliance operator in Voigt notation
        Cp = (C+abs(C))./2;
        Cm = (C-abs(C))./2;
        Dp = D'*Cp*D; % damaged part of stiffness operator in Voigt notation
        Dm = D'*Cm*D; % undamaged part of stiffness operator in Voigt notation
    else
        Dp = (D+abs(D))./2;
        Dm = (D-abs(D))./2;
    end
else
    if strcmpi(split,'stress')
        se = D*se; % stress tensor in Voigt notation
    end
    [Pp,Pm] = calc_proj_Miehe(mat,elem,xnode,xgauss,se); % projectors for strain tensor in Voigt notation
    
    P = calc_proj_notation(elem);
    P2 = P*P;
    if strcmpi(split,'stress')
        Pp = Pp*P2; % projector on positive part for stress tensor in Voigt notation
        Pm = Pm*P2; % projector on negative part for stress tensor in Voigt notation
        Dp = Pp*D; % damaged part of stiffness operator in Voigt notation
        Dm = Pm*D; % undamaged part of stiffness operator in Voigt notation
    else
        Pp = P2*Pp; % projector on positive part for strain tensor in Voigt notation
        Pm = P2*Pm; % projector on negative part for strain tensor in Voigt notation
        Dp = D*Pp; % damaged part of stiffness operator in Voigt notation
        Dm = D*Pm; % undamaged part of stiffness operator in Voigt notation
    end
end

%% Check stiffness/compliance operator decomposition
if ischarin('check',varargin)
    tol = 1e-12;
    D = calc_opmat(mat,elem,xnode,xgauss); % stiffness operator in Voigt notation
    if strcmpi(split,'stress')
        decompD = max(norm(D - (Dp+Dm))/norm(D),[],'all'); if decompD>tol, decompD, end
        decompDse = max(norm(D\se - (Dp+Dm)\se)/norm(D\se),[],'all'); if decompDse>tol, decompDse, end
        decompseDse = max(abs(se'/D*se - se'/(Dp+Dm)*se)/abs(se'/D*se),[],'all'); if decompseDse>tol, decompseDse, end
    else
        decompD = max(norm(D - (Dp+Dm))/norm(D),[],'all'); if decompD>tol, decompD, end
        decompDse = max(norm(D*se - (Dp+Dm)*se)/norm(D*se),[],'all'); if decompDse>tol, decompDse, end
        decompseDse = max(abs(se'*D*se - se'*(Dp+Dm)*se)/abs(se'*D*se),[],'all'); if decompseDse>tol, decompseDse, end
    end
end
