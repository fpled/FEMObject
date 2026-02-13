function [Dp,Dm] = calc_opmat_Miehe(mat,elem,xnode,xgauss,se,varargin)
% function [Dp,Dm] = calc_opmat_Miehe(mat,elem,xnode,xgauss,se)

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
    if strcmpi(split,'stress')
        Pp = Pp*P;
        Pm = Pm*P;
        if ischarin('+',varargin)
            % damaged positive part only
            Cp = Pp'/D*Pp; % damaged part of compliance operator in Kelvin-Mandel notation
            Cm = Pm'/D*Pm + Pp'/D*Pm + Pm'/D*Pp; % undamaged part of compliance operator in Kelvin-Mandel notation
        else
            % undamaged negative part only
            Cp = Pp'/D*Pp + Pp'/D*Pm + Pm'/D*Pp; % damaged part of compliance operator in Kelvin-Mandel notation
            Cm = Pm'/D*Pm; % undamaged part of compliance operator in Kelvin-Mandel notation
        end
        Cp = P'*Cp*P; % damaged part of compliance operator in Voigt notation
        Cm = P'*Cm*P; % undamaged part of compliance operator in Voigt notation
    else
        D = P'*D*P; % stiffness operator in Kelvin-Mandel notation
        Pp = P'*Pp;
        Pm = P'*Pm;
        if ischarin('+',varargin)
            % damaged positive part only
            Dp = Pp'*D*Pp; % damaged part of stiffness operator in Voigt notation
            Dm = Pm'*D*Pm + Pp'*D*Pm + Pm'*D*Pp; % undamaged part of stiffness operator in Voigt notation
        else
            % undamaged negative part only
            Dp = Pp'*D*Pp + Pp'*D*Pm + Pm'*D*Pp; % damaged part of stiffness operator in Voigt notation
            Dm = Pm'*D*Pm; % undamaged part of stiffness operator in Voigt notation
        end
    end
end

if strcmpi(split,'stress')
    Dp = D'*Cp*D; % damaged part of stiffness operator in Voigt notation
    Dm = D'*Cm*D; % undamaged part of stiffness operator in Voigt notation
end

%% Check stiffness/compliance operator decomposition
if ischarin('check',varargin)
    tol = 1e-12;
    D = calc_opmat(mat,elem,xnode,xgauss); % stiffness operator in Voigt notation
    if strcmpi(split,'stress')
        % C = inv(D); % compliance operator in Voigt notation
        % decompC = max(norm(C - (Cp+Cm))/norm(C),[],'all'); if decompC>tol, decompC, end
        % decompCse = max(norm(C*se - (Cp+Cm)*se)/norm(C*se),[],'all'); if decompCse>tol, decompCse, end
        % decompseCse = max(abs(se'*C*se - se'*(Cp+Cm)*se)/abs(se'*C*se),[],'all'); if decompseCse>tol, decompseCse, end
        decompD = max(norm(D - (Dp+Dm))/norm(D),[],'all'); if decompD>tol, decompD, end
        decompDse = max(norm(D\se - (Dp+Dm)\se)/norm(D\se),[],'all'); if decompDse>tol, decompDse, end
        decompseDse = max(abs(se'/D*se - se'/(Dp+Dm)*se)/abs(se'/D*se),[],'all'); if decompseDse>tol, decompseDse, end
    else
        decompD = max(norm(D - (Dp+Dm))/norm(D),[],'all'); if decompD>tol, decompD, end
        decompDse = max(norm(D*se - (Dp+Dm)*se)/norm(D*se),[],'all'); if decompDse>tol, decompDse, end
        decompseDse = max(abs(se'*D*se - se'*(Dp+Dm)*se)/abs(se'*D*se),[],'all'); if decompseDse>tol, decompseDse, end
    end
end
