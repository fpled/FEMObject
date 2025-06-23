function [Dp,Dm] = calc_opmat_He(mat,elem,xnode,xgauss,se,varargin)
% function [Dp,Dm] = calc_opmat_He(mat,elem,xnode,xgauss,se)

if nargin<=2
    xnode = [];
    xgauss = [];
    se = [];
end

dim = getdim(elem);
model = getparam(mat,'PFM'); % phase field model
split = getparam(mat,'PFS'); % phase field split

D = calc_opmat(mat,elem,xnode,xgauss); % stiffness operator in Voigt notation

P = calc_proj_notation(elem);
B = P*P';

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
        % C = inv(D); % compliance operator in Voigt notation
        % Cp = B'*(Pp'*C*Pp)*B; % damaged part of compliance operator in Voigt notation
        % Cm = B'*(Pm'*C*Pm)*B; % undamaged part of compliance operator in Voigt notation
        
        Cp = B'*(Pp'/D*Pp)*B; % damaged part of compliance operator in Voigt notation
        Cm = B'*(Pm'/D*Pm)*B; % undamaged part of compliance operator in Voigt notation
    else
        Dp = Pp'*(B'*D*B)*Pp; % damaged part of stiffness operator in Voigt notation
        Dm = Pm'*(B'*D*B)*Pm; % undamaged part of stiffness operator in Voigt notation
    end
end

if strcmpi(split,'stress')
    Dp = D'*Cp*D; % damaged part of stiffness operator in Voigt notation
    Dm = D'*Cm*D; % undamaged part of stiffness operator in Voigt notation
end

%% Check stiffness/compliance operator decomposition
% if verLessThan('matlab','9.1') % compatibility (<R2016b)
%     contain = @(str,pat) ~isempty(strfind(lower(str),pat));
% else
%     contain = @(str,pat) contains(str,pat,'IgnoreCase',true);
% end
% tol = 1e-12;
% D = calc_opmat(mat,elem,xnode,xgauss); % stiffness operator in Voigt notation
% if strcmpi(split,'stress')
%     % C = inv(D); % compliance operator in Voigt notation
%     % Cpp = D\(Pp*B); % damaged part of compliance operator in Voigt notation
%     % Cmm = D\(Pm*B); % undamaged part of compliance operator in Voigt notation
%     % decompC = max(norm(C - (Cpp+Cmm))/norm(C),[],'all'); if decompC>tol, decompC, end
%     % if contain(model,'amor')
%     %     decompC = max(norm(C - (Cp+Cm))/norm(C),[],'all'); if decompC>tol, decompC, end
%     % end
%     % decompCse = max(norm(C*se - (Cp+Cm)*se)/norm(C*se),[],'all'); if decompCse>tol, decompCse, end
%     % decompseCse = max(abs(se'*C*se - se'*(Cp+Cm)*se)/abs(se'*C*se),[],'all'); if decompseCse>tol, decompseCse, end
%     Dpp = Pp'*B'*D; % damaged part of stiffness operator in Voigt notation
%     Dmm = Pm'*B'*D; % undamaged part of stiffness operator in Voigt notation
%     decompD = max(norm(D - (Dpp+Dmm))/norm(D),[],'all'); if decompD>tol, decompD, end
%     decompDse = max(norm(D\se - (Dp+Dm)\se)/norm(D\se),[],'all'); if decompDse>tol, decompDse, end
%     decompseDse = max(abs(se'/D*se - se'/(Dp+Dm)*se)/abs(se'/D*se),[],'all'); if decompseDse>tol, decompseDse, end
% else
%     Dpp = D*B*Pp; % damaged part of stiffness operator in Voigt notation
%     Dmm = D*B*Pm; % undamaged part of stiffness operator in Voigt notation
%     decompD = max(norm(D - (Dpp+Dmm))/norm(D),[],'all'); if decompD>tol, decompD, end
%     if contain(model,'amor')
%         decompD = max(norm(D - (Dp+Dm))/norm(D),[],'all'); if decompD>tol, decompD, end
%     end
%     decompDse = max(norm(D*se - (Dp+Dm)*se)/norm(D*se),[],'all'); if decompDse>tol, decompDse, end
%     decompseDse = max(abs(se'*D*se - se'*(Dp+Dm)*se)/abs(se'*D*se),[],'all'); if decompseDse>tol, decompseDse, end
% end
