function [Dp,Dm] = calc_opmat_Miehe(mat,elem,xnode,xgauss,se,varargin)
% function [Dp,Dm] = calc_opmat_Miehe(mat,elem,xnode,xgauss,se)

if nargin<=2
    xnode = [];
    xgauss = [];
    se = [];
end

dim = getdim(elem);
split = getparam(mat,'PFS'); % phase field split

if dim>=2
    if strcmpi(split,'stress')
        D = calc_opmat(mat,elem,xnode,xgauss); % stiffness operator in Voigt notation
        se = D*se; % stress tensor in Voigt notation
    end
    [Pp,Pm] = calc_proj_Miehe(mat,elem,xnode,xgauss,se); % projectors for strain tensor in Voigt notation
    [Rp,Rm] = calc_proj_tr(mat,elem,xnode,xgauss,se); % projectors for trace of strain/stress tensor
    if strcmpi(split,'stress')
        P = calc_proj_notation(elem);
        Pp = P'*Pp*P; % projector on positive part in Kelvin-Mandel notation
        Pm = P'*Pm*P; % projector on negative part in Kelvin-Mandel notation
    end
end

switch dim
    case 1
        D = calc_opmat(mat,elem,xnode,xgauss); % stiffness operator in Voigt notation
        if strcmpi(split,'stress')
            C = inv(D); % compliance operator in Voigt notation
            Cp = (C+abs(C))./2;
            Cm = (C-abs(C))./2;
        else
            Dp = (D+abs(D))./2;
            Dm = (D-abs(D))./2;
        end
    case 2
        E = evalparam(mat,'E',elem,xnode,xgauss); % Young modulus
        nu = evalparam(mat,'NU',elem,xnode,xgauss); % Poisson ratio
        mu = E/(1+nu)/2; % second Lame coefficient (shear modulus)
        if isaxi(elem)
            v = [1,1,1,0]';
            I = v*v';
            if strcmpi(split,'stress')
                Cp = -nu/E*Rp*I + 1/(2*mu)*Pp; % damaged part of compliance operator in Kelvin-Mandel notation
                Cm = -nu/E*Rm*I + 1/(2*mu)*Pm; % undamaged part of compliance operator in Kelvin-Mandel notation
                Cp = P'*Cp*P; % damaged part of compliance operator in Voigt notation
                Cm = P'*Cm*P; % undamaged part of compliance operator in Voigt notation
            else
                lambda = E*nu/(1+nu)/(1-2*nu); % first Lame coefficient
                Dp = lambda*Rp*I + 2*mu*Pp; % damaged part of stiffness operator in Voigt notation
                Dm = lambda*Rm*I + 2*mu*Pm; % undamaged part of stiffness operator in Voigt notation
            end
        else
            e = evalparam(mat,'DIM3',elem,xnode,xgauss); % thickness
            v = [1,1,0]';
            I = v*v';
            if strcmpi(split,'stress')
                if strcmp(getoption(elem),'DEFO')
                    nu = nu*(1+nu);
                end
                Cp = 1/e*(-nu/E*Rp*I + 1/(2*mu)*Pp); % damaged part of compliance operator in Kelvin-Mandel notation
                Cm = 1/e*(-nu/E*Rm*I + 1/(2*mu)*Pm); % undamaged part of compliance operator in Kelvin-Mandel notation
                Cp = P'*Cp*P; % damaged part of compliance operator in Voigt notation
                Cm = P'*Cm*P; % undamaged part of compliance operator in Voigt notation
            else
                switch getoption(elem)
                    case 'DEFO'
                        lambda = E*nu/(1+nu)/(1-2*nu); % first Lame coefficient
                    otherwise
                        lambda = E*nu/(1-nu^2); % first Lame coefficient
                end
                Dp = e*(lambda*Rp*I + 2*mu*Pp); % damaged part of stiffness operator in Voigt notation
                Dm = e*(lambda*Rm*I + 2*mu*Pm); % undamaged part of stiffness operator in Voigt notation
            end
        end
        
    case 3
        E = evalparam(mat,'E',elem,xnode,xgauss); % Young modulus
        nu = evalparam(mat,'NU',elem,xnode,xgauss); % Poisson ratio
        mu = E/(1+nu)/2; % second Lame coefficient (shear modulus)
        v = [1,1,1,0,0,0]';
        I = v*v';
        if strcmpi(split,'stress')
            Cp = -nu/E*Rp*I + 1/(2*mu)*Pp; % damaged part of compliance operator in Kelvin-Mandel notation
            Cm = -nu/E*Rm*I + 1/(2*mu)*Pm; % undamaged part of compliance operator in Kelvin-Mandel notation
            Cp = P'*Cp*P; % damaged part of compliance operator in Voigt notation
            Cm = P'*Cm*P; % undamaged part of compliance operator in Voigt notation
        else
            lambda = E*nu/(1+nu)/(1-2*nu); % first Lame coefficient
            Dp = lambda*Rp*I + 2*mu*Pp; % damaged part of stiffness operator in Voigt notation
            Dm = lambda*Rm*I + 2*mu*Pm; % undamaged part of stiffness operator in Voigt notation
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
