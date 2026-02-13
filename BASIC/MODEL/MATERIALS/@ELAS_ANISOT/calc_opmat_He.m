function [Dp,Dm] = calc_opmat_He(mat,elem,xnode,xgauss,se,varargin)
% function [Dp,Dm] = calc_opmat_He(mat,elem,xnode,xgauss,se)

if nargin<=2
    xnode = [];
    xgauss = [];
    se = [];
end

dim = getdim(elem);
split = getparam(mat,'PFS'); % phase field split

D = calc_opmat(mat,elem,xnode,xgauss); % stiffness operator in Voigt notation

P = calc_proj_notation(elem);
P2 = P*P;

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
    [Pp,Pm] = calc_proj_He(mat,elem,xnode,xgauss,se); % projectors for strain tensor in Voigt notation
    
    if strcmpi(split,'stress')
        Pp = Pp*P2; % projector on positive part for stress tensor in Voigt notation
        Pm = Pm*P2; % projector on negative part for stress tensor in Voigt notation
        % C = inv(D); % compliance operator in Voigt notation
        % % Cp = Pp'*C*Pp; % secant representation of damaged part of compliance operator in Voigt notation
        % % Cm = Pm'*C*Pm; % secant representation of undamaged part of compliance operator in Voigt notation
        % Cp = C*Pp; % tangent representation of damaged part of compliance operator in Voigt notation
        % Cm = C*Pm; % tangent representation of undamaged part of compliance operator in Voigt notation
        
        % % Cp = Pp'/D*Pp; % secant representation of damaged part of compliance operator in Voigt notation
        % % Cm = Pm'/D*Pm; % secant representation of undamaged part of compliance operator in Voigt notation
        % Cp = D\Pp; % tangent representation of damaged part of compliance operator in Voigt notation
        % Cm = D\Pm; % tangent representation of undamaged part of compliance operator in Voigt notation
        
        % Dp = D'*Cp*D; % damaged part of stiffness operator in Voigt notation
        % Dm = D'*Cm*D; % undamaged part of stiffness operator in Voigt notation

        Dp = Pp*D; % damaged part of stiffness operator in Voigt notation
        Dm = Pm*D; % undamaged part of stiffness operator in Voigt notation
    else
        Pp = P2*Pp; % projector on positive part for strain tensor in Voigt notation
        Pm = P2*Pm; % projector on negative part for strain tensor in Voigt notation
        % Dp = Pp'*D*Pp; % secant representation of damaged part of stiffness operator in Voigt notation
        % Dm = Pm'*D*Pm; % secant representation of undamaged part of stiffness operator in Voigt notation
        Dp = D*Pp; % tangent representation of damaged part of stiffness operator in Voigt notation
        Dm = D*Pm; % tangent representation of undamaged part of stiffness operator in Voigt notation
    end
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
