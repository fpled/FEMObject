function [Dp,Dm] = calc_opmat_Amor(mat,elem,xnode,xgauss,se,varargin)
% function [Dp,Dm] = calc_opmat_Amor(mat,elem,xnode,xgauss,se)

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
    [Rp,Rm] = calc_proj_tr(mat,elem,xnode,xgauss,se); % projectors for trace of strain/stress tensor
end

switch dim
    case 1
        D = calc_opmat(mat,elem,xnode,xgauss);
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
        mu = E/(1+nu)/2; % second Lamé coefficient (shear modulus)
        if isaxi(elem)
            lambda = E*nu/(1+nu)/(1-2*nu); % first Lamé coefficient
            kappa = lambda+mu; % bulk modulus
            v = [1,1,1,0]';
            I = v*v';
            if strcmpi(split,'stress')
                Id = diag([1,1,1,2]);
                Pd = Id-I/2;
                Cp = 1/(9*kappa)*Rp*I + 1/(2*mu)*Pd; % damaged part of compliance operator in Voigt notation
                Cm = 1/(9*kappa)*Rm*I; % undamaged part of compliance operator in Voigt notation
            else
                Id = diag([1,1,1,1/2]);
                Pd = Id-I/2;
                Dp = kappa*Rp*I + 2*mu*Pd; % damaged part of stiffness operator in Voigt notation
                Dm = kappa*Rm*I; % undamaged part of stiffness operator in Voigt notation
            end
        else
            e = evalparam(mat,'DIM3',elem,xnode,xgauss); % thickness
            switch getoption(elem)
                case 'DEFO'
                    lambda = E*nu/(1+nu)/(1-2*nu); % first Lamé coefficient
                otherwise
                    lambda = E*nu/(1-nu^2); % first Lamé coefficient
            end
            kappa = lambda+mu; % bulk modulus
            v = [1,1,0]';
            I = v*v';
            if strcmpi(split,'stress')
                Id = diag([1,1,2]);
                Pd = Id-I/2;
                Cp = 1/e*(1/(4*kappa)*Rp*I + 1/(2*mu)*Pd); % damaged part of compliance operator in Voigt notation
                Cm = 1/e*(1/(4*kappa)*Rm*I); % undamaged part of compliance operator in Voigt notation
            else
                Id = diag([1,1,1/2]);
                Pd = Id-I/2;
                Dp = e*(kappa*Rp*I + 2*mu*Pd); % damaged part of stiffness operator in Voigt notation
                Dm = e*(kappa*Rm*I); % undamaged part of stiffness operator in Voigt notation
            end
        end
        
    case 3
        E = evalparam(mat,'E',elem,xnode,xgauss); % Young modulus
        nu = evalparam(mat,'NU',elem,xnode,xgauss); % Poisson ratio
        lambda = E*nu/(1+nu)/(1-2*nu); % first Lamé coefficient
        mu = E/(1+nu)/2; % second Lamé coefficient (shear modulus)
        kappa = lambda+2*mu/3; % bulk modulus
        v = [1,1,1,0,0,0]';
        I = v*v';
        if strcmpi(split,'stress')
            Id = diag([1,1,1,2,2,2]);
            Pd = Id-I/3;
            Cp = 1/(9*kappa)*Rp*I + 1/(2*mu)*Pd; % damaged part of compliance operator in Voigt notation
            Cm = 1/(9*kappa)*Rm*I; % undamaged part of compliance operator in Voigt notation
        else
            Id = diag([1,1,1,1/2,1/2,1/2]);
            Pd = Id-I/3;
            Dp = kappa*Rp*I + 2*mu*Pd; % damaged part of stiffness operator in Voigt notation
            Dm = kappa*Rm*I; % undamaged part of stiffness operator in Voigt notation
        end
        
end

% Other implementation
if dim>=2
    [Pp,Pm] = calc_proj_Amor(mat,elem,xnode,xgauss,se); % projectors for strain tensor in Voigt notation
    
    P = calc_proj_notation(elem);
    P2 = P*P;
    if strcmpi(split,'stress')
        Pp = Pp*P2; % projector on positive part for stress tensor in Voigt notation
        Pm = Pm*P2; % projector on negative part for stress tensor in Voigt notation
        Cp = Pp'/D*Pp; % damaged part of compliance operator in Voigt notation
        Cm = Pm'/D*Pm; % undamaged part of compliance operator in Voigt notation
    else
        Pp = P2*Pp; % projector on positive part for strain tensor in Voigt notation
        Pm = P2*Pm; % projector on negative part for strain tensor in Voigt notation
        D = calc_opmat(mat,elem,xnode,xgauss); % stiffness operator in Voigt notation
        Dp = Pp'*D*Pp; % damaged part of stiffness operator in Voigt notation
        Dm = Pm'*D*Pm; % undamaged part of stiffness operator in Voigt notation
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
