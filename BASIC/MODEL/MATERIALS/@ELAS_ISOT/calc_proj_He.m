function [Pp,Pm] = calc_proj_He(mat,elem,xnode,xgauss,se,varargin)
% function [Pp,Pm] = calc_proj_He(mat,elem,xnode,xgauss,se)

model = getparam(mat,'PFM'); % phase field model
split = getparam(mat,'PFS'); % phase field split

P = calc_proj_notation(elem);
if strcmpi(split,'stress')
    B = P*P';
    se = B*se; % stress tensor converted as strain tensor in Voigt notation
end

%% Numerical computation of the square root of stiffness/compliance tensor
D = calc_opmat(mat,elem,xnode,xgauss); % stiffness operator in Voigt notation
D = P'*D*P; % stiffness operator in Kelvin-Mandel notation
sqrtD = sqrtm(D); % square root of stiffness operator in Kelvin-Mandel notation
% sqrtC = inv(sqrtD); % square root of compliance operator in Kelvin-Mandel notation

%% Analytical computation of the square root of stiffness/compliance tensor
% [sqrtDana,sqrtCana] = calc_sqrtopmat(mat,elem,xnode,xgauss); % square root of stiffness operator in Voigt notation
% sqrtDana = P*sqrtDana; % square root of stiffness operator in Kelvin-Mandel notation
% sqrtCana = P\sqrtCana; % square root of compliance operator in Kelvin-Mandel notation

%% Check numerical and analytical computations
% tol = 1e-12;
% errsqrtD = max(norm(sqrtDana-sqrtD)/norm(sqrtDana),[],'all'); if errsqrtD>tol, keyboard, end
% errsqrtC = max(norm(sqrtCana-sqrtC)/norm(sqrtCana),[],'all'); if errsqrtC>tol, keyboard, end

%% Transformed strain/stress tensor
if strcmpi(split,'stress')
    % sqrtC = inv(sqrtD); % square root of compliance operator in Kelvin-Mandel notation
    % set = P\sqrtC*(P\se);
    % set = (P\sqrtC/P)*se;
    set = P\(sqrtD\(P\se)); % transformed stress tensor in Voigt notation
else
    % set = P*sqrtD*(P\se);
    set = (P*sqrtD/P)*se; % transformed strain tensor in Voigt notation
end

%% Computation of projectors of strain/stress tensor onto volumetric/spherical/hydrostatic expansion (or tensile bulk) and deviatoric parts for generalized Amor decomposition
%% Computation of projectors of strain/stress tensor onto positive (tensile) and negative (compressive) parts for generalized Freddi decomposition
% Projectors for transformed strain/stress tensor
if verLessThan('matlab','9.1') % compatibility (<R2016b)
    projAmor = ~isempty(strfind(lower(model),'amor'));
    projFreddi = ~isempty(strfind(lower(model),'freddi'));
else
    projAmor = contains(model,'amor','IgnoreCase',true);
    projFreddi = contains(model,'freddi','IgnoreCase',true);
end
if projAmor
    [Ppt,Pmt] = calc_proj_Amor(mat,elem,xnode,xgauss,set); % projectors for strain tensor in Voigt notation
elseif projFreddi
    [Ppt,Pmt] = calc_proj_Miehe(mat,elem,xnode,xgauss,set); % projectors for strain tensor in Voigt notation
else
    error(['Wrong phase field model ' model])
end
Ppt = P'*Ppt*P; % projector on positive part in Kelvin-Mandel notation
Pmt = P'*Pmt*P; % projector on negative part in Kelvin-Mandel notation

% Projectors for strain/stress tensor
if strcmpi(split,'stress')
    % Pp = sqrtC\(Ppt*sqrtC);
    % Pm = sqrtC\(Pmt*sqrtC);
    % Pp = sqrtC\Ppt*sqrtC;
    % Pm = sqrtC\Pmt*sqrtC;
    Pp = sqrtD*Ppt/sqrtD; % projector on positive part in Kelvin-Mandel notation
    Pm = sqrtD*Pmt/sqrtD; % projector on negative part in Kelvin-Mandel notation
else
    Pp = sqrtD\(Ppt*sqrtD); % projector on positive part in Kelvin-Mandel notation
    Pm = sqrtD\(Pmt*sqrtD); % projector on negative part in Kelvin-Mandel notation
end
Pp = P'\Pp/P; % projector on positive part in Voigt notation
Pm = P'\Pm/P; % projector on negative part in Voigt notation

%% Check orthogonality condition and strain/stress decomposition
% tol = 1e-12;
% sep = P*(Pp*se); % positive part of strain/stress tensor in Kelvin-Mandel notation
% sem = P*(Pm*se); % negative part of strain/stress tensor in Kelvin-Mandel notation
% se = P\se; % strain/stress tensor in Kelvin-Mandel notation
% if strcmpi(split,'stress')
%     % C = inv(D); % compliance operator in Kelvin-Mandel notation
%     % orthCpm = max(abs(sep'*C*sem)/abs(se'*C*se),[],'all'); if orthCpm>tol, orthCpm, end
%     % orthCmp = max(abs(sem'*C*sep)/abs(se'*C*se),[],'all'); if orthCmp>tol, orthCmp, end
%     orthDpm = max(abs(sep'/D*sem)/abs(se'/D*se),[],'all'); if orthDpm>tol, orthDpm, end
%     orthDmp = max(abs(sem'/D*sep)/abs(se'/D*se),[],'all'); if orthDmp>tol, orthDmp, end
% else
%     orthDpm = max(abs(sep'*D*sem)/abs(se'*D*se),[],'all'); if orthDpm>tol, orthDpm, end
%     orthDmp = max(abs(sem'*D*sep)/abs(se'*D*se),[],'all'); if orthDmp>tol, orthDmp, end
% end
% decomp = max(norm(se - (sep+sem))/norm(se),[],'all'); if decomp>tol, decomp, end

end
