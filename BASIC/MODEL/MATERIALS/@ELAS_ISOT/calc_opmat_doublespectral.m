function [Dp,Dm] = calc_opmat_doublespectral(mat,elem,xnode,xgauss,varargin)
% function [Dp,Dm] = calc_opmat_doublespectral(mat,elem,xnode,xgauss)

if nargin<=2
    xnode = [];
    xgauss = [];
end

dim = getdim(elem);
% model = getparam(mat,'PFM'); % phase field model
% split = getparam(mat,'PFS'); % phase field split

D = calc_opmat(mat,elem,xnode,xgauss); % stiffness operator in Voigt notation

% P = calc_proj_notation(elem);
% P2 = P*P;

%% Numerical spectral (eigenvalue) decomposition of elasticity tensor
% Eigenvalues and eigenvectors
[Venum,Denum] = eig(MYDOUBLEND(D),'vector','sort');

% Eigenprojectors
% switch dim
%     case 1
%         M1num = Venum(:,1)*Venum(:,1)';
%     case 2
%         M1num = Venum(:,1)*Venum(:,1)';
%         M2num = Venum(:,2)*Venum(:,2)';
%         M3num = Venum(:,3)*Venum(:,3)';
%         if isaxi(elem)
%             M4num = Venum(:,4)*Venum(:,4)';
%         end
%     case 3
%         M1num = Venum(:,1)*Venum(:,1)';
%         M2num = Venum(:,2)*Venum(:,2)';
%         M3num = Venum(:,3)*Venum(:,3)';
%         M4num = Venum(:,4)*Venum(:,4)';
%         M5num = Venum(:,5)*Venum(:,5)';
%         M6num = Venum(:,6)*Venum(:,6)';
% end

%% Check numerical spectral (eigenvalue) decomposition of elasticity tensor
% tol = 1e-12;
% if verLessThan('matlab','9.1') % compatibility (<R2016b)
%     Dnum = Venum*diag(Denum)*Venum';
% else
%     Dnum = Venum*(Denum.*Venum');
% end
% decompDnum = max(norm(D - Dnum)/norm(D),[],'all'); if decompDnum>tol, decompDnum, end

%% Analytical spectral (eigenvalue) decomposition of elasticity tensor
% switch dim
%     case 1
%         De = D;
%         Ve = 1;
%         M1 = 1;
%     case 2
%         E = evalparam(mat,'E',elem,xnode,xgauss); % Young modulus
%         nu = evalparam(mat,'NU',elem,xnode,xgauss); % Poisson ratio
%         mu = E/(1+nu)/2; % second Lamé coefficient (shear modulus)
%         if isaxi(elem)
%             lambda = E*nu/(1+nu)/(1-2*nu); % first Lamé coefficient
%             kappa = E/(1-2*nu)/3; % bulk modulus
%             % kappa = lambda+2/3*mu; % bulk modulus
%             
%             % Eigenvalues (in ascending order)
%             De = [mu, 2*mu, 2*mu, 3*kappa]';
%             
%             % Eigenvectors
%             v1 = [0,0,0,1]';
%             v2 = 1/sqrt(5)*[sqrt(2)/2-1,sqrt(2)/2+1,-sqrt(2),0]';
%             v3 = 1/sqrt(15)*[1+3*sqrt(2)/2,1-3*sqrt(2)/2,-2,0]';
%             v4 = 1/sqrt(3)*[1,1,1,0]';
%             Ve = [v1, v2, v3, v4];
%             
%             % Eigenprojectors
%             M1 = v1*v1';
%             M2 = v2*v2';
%             M3 = v3*v3';
%             M4 = v4*v4';
%             % M1 = [0,0,0,0; 0,0,0,0; 0,0,0,0; 0,0,0,1];
%             % M2 = 1/5*[[3/2-sqrt(2),-1/2,-(1-sqrt(2)); -1/2,3/2+sqrt(2),-(1+sqrt(2)); -(1-sqrt(2)),-(1+sqrt(2)),2], zeros(3,1); zeros(1,4)];
%             % M3 = 1/15*[[11/2+3*sqrt(2),-7/2,-(2+3*sqrt(2)); -7/2,11/2-3*sqrt(2),-(2-3*sqrt(2)); -(2+3*sqrt(2)),-(2-3*sqrt(2)),4], zeros(3,1); zeros(1,4)];
%             % M4 = 1/3*[ones(3), zeros(3,1); zeros(1,4)];
%         else
%             e = evalparam(mat,'DIM3',elem,xnode,xgauss); % thickness
%             switch getoption(elem)
%                 case 'DEFO'
%                     lambda = E*nu/(1+nu)/(1-2*nu); % first Lamé coefficient
%                     % kappa = E/(1+nu)/(1-2*nu)/2; % bulk modulus
%                 otherwise
%                     lambda = E*nu/(1-nu^2); % first Lamé coefficient
%                     % kappa = E/(1-nu)/2; % bulk modulus
%             end
%             kappa = lambda+mu; % bulk modulus
%             
%             % Eigenvalues (in ascending order)
%             De = e*[mu, 2*mu, 2*kappa]';
%             
%             % Eigenvectors
%             v1 = [0,0,1]';
%             v2 = 1/sqrt(2)*[-1,1,0]';
%             v3 = 1/sqrt(2)*[1,1,0]';
%             Ve = [v1, v2, v3];
%             
%             % Eigenprojectors
%             M1 = v1*v1';
%             M2 = v2*v2';
%             M3 = v3*v3';
%             % M1 = [0,0,0; 0,0,0; 0,0,1];
%             % M2 = 1/2*[1,-1,0; -1,1,0; 0,0,0];
%             % M3 = 1/2*[1,1,0; 1,1,0; 0,0,0];
%         end
%         
%     case 3
%         E = evalparam(mat,'E',elem,xnode,xgauss); % Young modulus
%         nu = evalparam(mat,'NU',elem,xnode,xgauss); % Poisson ratio
%         lambda = E*nu/(1+nu)/(1-2*nu); % first Lamé coefficient
%         mu = E/(1+nu)/2; % second Lamé coefficient (shear modulus)
%         kappa = E/(1-2*nu)/3; % bulk modulus
%         % kappa = lambda+2/3*mu; % bulk modulus
%         
%         % Eigenvalues (in ascending order)
%         De = [mu*ones(1,3), 2*mu*ones(1,2), 3*kappa]';
%         
%         % Eigenvectors
%         v1 = [0,0,0,1,0,0]';
%         v2 = [0,0,0,0,1,0]';
%         v3 = [0,0,0,0,0,1]';
%         v4 = 1/sqrt(5)*[1/sqrt(2)-1,1/sqrt(2)+1,-sqrt(2),0,0,0]';
%         v5 = 1/sqrt(15)*[1+3/sqrt(2),1-3/sqrt(2),-2,0,0,0]';
%         v6 = 1/sqrt(3)*[1,1,1,0,0,0]';
%         Ve = [v1, v2, v3, v4, v5, v6];
%         
%         % Eigenprojectors
%         M1 = v1*v1';
%         M2 = v2*v2';
%         M3 = v3*v3';
%         M4 = v4*v4';
%         M5 = v5*v5';
%         M6 = v6*v6';
%         % M1 = [zeros(3,6); zeros(3), [1,0,0; 0,0,0; 0,0,0]];
%         % M2 = [zeros(3,6); zeros(3), [0,0,0; 0,1,0; 0,0,0]];
%         % M3 = [zeros(3,6); zeros(3), [0,0,0; 0,0,0; 0,0,1]];
%         % M4 = 1/5*[[3/2-sqrt(2),-1/2,-(1-sqrt(2)); -1/2,3/2+sqrt(2),-(1+sqrt(2)); -(1-sqrt(2)),-(1+sqrt(2)),2], zeros(3); zeros(3,6)];
%         % M5 = 1/15*[[11/2+3*sqrt(2),-7/2,-(2+3*sqrt(2)); -7/2,11/2-3*sqrt(2),-(2-3*sqrt(2)); -(2+3*sqrt(2)),-(2-3*sqrt(2)),4], zeros(3); zeros(3,6)];
%         % M6 = 1/3*[ones(3), zeros(3); zeros(3,6)];
%         
% end

%% Check analytical spectral (eigenvalue) decomposition of elasticity tensor
% tol = 1e-12;
% if verLessThan('matlab','9.1') % compatibility (<R2016b)
%     Dana = Ve*diag(De)*Ve';
% else
%     Dana = Ve*(De.*Ve');
% end
% decompDana = max(norm(D - Dana)/norm(D),[],'all'); if decompDana>tol, decompDana, end

%% Check numerical and analytical spectral (eigenvalue) decompositions of elasticity tensor
% tol = 1e-12;
% switch dim
%     case 1
%         Dnum = Denum(1)*M1num;
%         Dana = De(1)*M1;
%     case 2
%         if isaxi(elem)
%             Dnum = Denum(1)*M1num + Denum(2)*M2num + Denum(3)*M3num + Denum(4)*M4num;
%             Dana = De(1)*M1 + De(2)*M2 + De(3)*M3 + De(4)*M4;
%         else
%             Dnum = Denum(1)*M1num + Denum(2)*M2num + Denum(3)*M3num;
%             Dana = De(1)*M1 + De(2)*M2 + De(3)*M3;
%         end
%     case 3
%         Dnum = Denum(1)*M1num + Denum(2)*M2num + Denum(3)*M3num...
%             + Denum(4)*M4num + Denum(5)*M5num + Denum(6)*M6num;
%         Dana = De(1)*M1 + De(2)*M2 + De(3)*M3...
%             + De(4)*M4 + De(5)*M5 + De(6)*M6;
% end
% errnum = max(norm(D-Dnum)/norm(D),[],'all'); if errnum>tol, errnum, end
% errana = max(norm(D-Dana)/norm(D),[],'all'); if errana>tol, errana, end
% 
% errval = max(abs(De-Denum)./abs(De),[],'all'); if errval>tol, errval, end
% errM1 = max(norm(M1-M1num)/norm(M1),[],'all'); if errM1>tol, errM1, end
% if dim>=2
%     errM2 = max(norm(M2-M2num)/norm(M2),[],'all'); if errM2>tol, errM2, end
%     errM3 = max(norm(M3-M3num)/norm(M3),[],'all'); if errM3>tol, errM3, end
%     if isaxi(elem)
%         errM4 = max(norm(M4-M4num)/norm(M4),[],'all'); if errM4>tol, errM4, end
%     end
% end
% if dim==3
%     errM4 = max(norm(M4-M4num)/norm(M4),[],'all'); if errM4>tol, errM4, end
%     errM5 = max(norm(M5-M5num)/norm(M5),[],'all'); if errM5>tol, errM5, end
%     errM6 = max(norm(M6-M6num)/norm(M6),[],'all'); if errM6>tol, errM6, end
% end
% switch dim
%     case 1
%         errM = max(norm(M1-eye(1)),[],'all'); if errM>tol, errM, end
%         errMnum = max(norm(M1num-eye(1)),[],'all'); if errMnum>tol, errMnum, end
%     case 2
%         if isaxi(elem)
%             errM = max(norm(M1+M2+M3+M4-eye(4)),[],'all'); if errM>tol, errM, end
%             errMnum = max(norm(M1num+M2num+M3num+M4num-eye(4)),[],'all'); if errMnum>tol, errMnum, end
%         else
%             errM = max(norm(M1+M2+M3-eye(3)),[],'all'); if errM>tol, errM, end
%             errMnum = max(norm(M1num+M2num+M3num-eye(3)),[],'all'); if errMnum>tol, errMnum, end
%         end
%     case 3
%         errM = max(norm(M1+M2+M3+M4+M5+M6-eye(6)),[],'all'); if errM>tol, errM, end
%         errMnum = max(norm(M1num+M2num+M3num+M4num+M5num+M6num-eye(6)),[],'all'); if errMnum>tol, errMnum, end
% end

%% Numerical spectral (eigenvalue) decomposition of eigenvectors
Ve = Venum;
De = Denum;
% Eigenvalues and eigenvectors
if isaxi(elem)
    L = zerosND([3,size(Ve,2),sizeND(Ve)]);
else
    L = zerosND([dim,size(Ve,2),sizeND(Ve)]);
end
L2 = zerosND([dim,dim,size(Ve,2),sizeND(Ve)]);
N = zerosND([size(Ve,1),dim,size(Ve,2),sizeND(Ve)]);
for i=1:size(Ve,2)
    % Eigenvector #i
    vi = Ve(:,i);
    % Eigenmatrix #i
    switch dim
        case 1
            Vi = vi(1);
        case 2
            if isaxi(elem)
                Vi = [vi(1) vi(4) 0;
                    vi(4) vi(2) 0;
                    0 0 vi(3)];
            else
                Vi = [vi(1) vi(3);
                    vi(3) vi(2)];
            end
        case 3
            Vi = [vi(1) vi(4) vi(5);
                vi(4) vi(2) vi(6);
                vi(5) vi(6) vi(3)];
    end
    
    % Eigenvalues and eigenvectors for eigenmatrix #i
    [Ui,Li] = eig(MYDOUBLEND(Vi),'vector','sort');
    
    %% Check numerical spectral (eigenvalue) decomposition of eigenmatrix #i
    % tol = 1e-12;
    % if verLessThan('matlab','9.1') % compatibility (<R2016b)
    %     Vinum = Ui*diag(Li)*Ui';
    % else
    %     Vinum = Ui*(Li.*Ui');
    % end
    % decompVi = max(norm(Vi - Vinum)/norm(Vi),[],'all'); if decompVi>tol, decompVi, end
    
    % Eigentensors for eigenvector #i
    ni = zerosND([size(vi,1),size(Ui,2),sizeND(Ui)]);
    for j=1:size(Ui,2)
        Nij = Ui(:,j)*Ui(:,j)';
        switch dim
            case 1
                nij = Nij(1,1);
            case 2
                if isaxi(elem)
                    nij = [Nij(1,1) Nij(2,2) Nij(3,3) Nij(1,2)]';
                else
                    nij = [Nij(1,1) Nij(2,2) Nij(1,2)]';
                end
            case 3
                nij = [Nij(1,1) Nij(2,2) Nij(3,3) Nij(1,2) Nij(1,3) Nij(2,3)]';
        end
        ni(:,j) = nij;
    end
    
    %% Check numerical spectral (eigenvalue) decomposition of eigenvector #i
    % tol = 1e-12;
    % switch dim
    %     case 1
    %         vinum = Li(1)*ni(:,1);
    %     case 2
    %         if isaxi(elem)
    %             vinum = Li(1)*ni(:,1) + Li(2)*ni(:,2) + Li(3)*ni(:,3);
    %         else
    %             vinum = Li(1)*ni(:,1) + Li(2)*ni(:,2);
    %         end
    %     case 3
    %         vinum = Li(1)*ni(:,1) + Li(2)*ni(:,2) + Li(3)*ni(:,3);
    % end
    % vinum = ni*Li; 
    % % vinum = sum(Li'.*ni);
    % decompvi = max(norm(vinum-vi)/norm(vi),[],'all'); if decompvi>tol, decompvi, end
    
    %% Check numerical spectral (eigenvalue) decomposition of eigentensor #i
    % Mi = vi*vi';
    % L2i = Li*Li';
    % Minum = ni*L2i*ni';
    % decompMi = max(norm(Minum-Mi)/norm(Mi),[],'all'); if decompMi>tol, decompMi, end
    
    N(:,:,i) = ni;
    L(:,i) = Li;
    L2(:,:,i) = Li*Li';
end

%% Check double spectral (eigenvalue) decomposition of elasticity tensor
% tol = 1e-12;
% switch dim
%     case 1
%         Ddouble = De(1)*N(:,:,1)*L2(:,:,1)*N(:,:,1)';
%     case 2
%         if isaxi(elem)
%             Ddouble = De(1)*N(:,:,1)*L2(:,:,1)*N(:,:,1)'...
%                 + De(2)*N(:,:,2)*L2(:,:,2)*N(:,:,2)'...
%                 + De(3)*N(:,:,3)*L2(:,:,3)*N(:,:,3)'...
%                 + De(4)*N(:,:,4)*L2(:,:,4)*N(:,:,4)';
%         else
%             Ddouble = De(1)*N(:,:,1)*L2(:,:,1)*N(:,:,1)'...
%                 + De(2)*N(:,:,2)*L2(:,:,2)*N(:,:,2)'...
%                 + De(3)*N(:,:,3)*L2(:,:,3)*N(:,:,3)';
%         end
%     case 3
%         Ddouble = De(1)*N(:,:,1)*L2(:,:,1)*N(:,:,1)'...
%             + De(2)*N(:,:,2)*L2(:,:,2)*N(:,:,2)'...
%             + De(3)*N(:,:,3)*L2(:,:,3)*N(:,:,3)'...
%             + De(4)*N(:,:,4)*L2(:,:,4)*N(:,:,4)'...
%             + De(5)*N(:,:,5)*L2(:,:,5)*N(:,:,5)'...
%             + De(6)*N(:,:,6)*L2(:,:,6)*N(:,:,6)';
% end
% decompD = max(norm(D - Ddouble)/norm(D),[],'all'); if decompD>tol, decompD, end

%% Double spectral decomposition of elasticity tensor
L2p = (L2+abs(L2))./2;
L2m = (L2-abs(L2))./2;
switch dim
    case 1
        Dp = De(1)*N(:,:,1)*L2p(:,:,1)*N(:,:,1)'; % damaged part of stiffness operator in Voigt notation
        Dm = De(1)*N(:,:,1)*L2m(:,:,1)*N(:,:,1)'; % undamaged part of stiffness operator in Voigt notation
    case 2
        if isaxi(elem)
            Dp = De(1)*N(:,:,1)*L2p(:,:,1)*N(:,:,1)'...
                + De(2)*N(:,:,2)*L2p(:,:,2)*N(:,:,2)'...
                + De(3)*N(:,:,3)*L2p(:,:,3)*N(:,:,3)'...
                + De(4)*N(:,:,4)*L2p(:,:,4)*N(:,:,4)'; % damaged part of stiffness operator in Voigt notation
            Dm = De(1)*N(:,:,1)*L2m(:,:,1)*N(:,:,1)'...
                + De(2)*N(:,:,2)*L2m(:,:,2)*N(:,:,2)'...
                + De(3)*N(:,:,3)*L2m(:,:,3)*N(:,:,3)'...
                + De(4)*N(:,:,4)*L2m(:,:,4)*N(:,:,4)'; % undamaged part of stiffness operator in Voigt notation
        else
            Dp = De(1)*N(:,:,1)*L2p(:,:,1)*N(:,:,1)'...
                + De(2)*N(:,:,2)*L2p(:,:,2)*N(:,:,2)'...
                + De(3)*N(:,:,3)*L2p(:,:,3)*N(:,:,3)'; % damaged part of stiffness operator in Voigt notation
            Dm = De(1)*N(:,:,1)*L2m(:,:,1)*N(:,:,1)'...
                + De(2)*N(:,:,2)*L2m(:,:,2)*N(:,:,2)'...
                + De(3)*N(:,:,3)*L2m(:,:,3)*N(:,:,3)'; % undamaged part of stiffness operator in Voigt notation
        end
    case 3
        Dp = De(1)*N(:,:,1)*L2p(:,:,1)*N(:,:,1)'...
            + De(2)*N(:,:,2)*L2p(:,:,2)*N(:,:,2)'...
            + De(3)*N(:,:,3)*L2p(:,:,3)*N(:,:,3)'...
            + De(4)*N(:,:,4)*L2p(:,:,4)*N(:,:,4)'...
            + De(5)*N(:,:,5)*L2p(:,:,5)*N(:,:,5)'...
            + De(6)*N(:,:,6)*L2p(:,:,6)*N(:,:,6)'; % damaged part of stiffness operator in Voigt notation
        Dm = De(1)*N(:,:,1)*L2m(:,:,1)*N(:,:,1)'...
            + De(2)*N(:,:,2)*L2m(:,:,2)*N(:,:,2)'...
            + De(3)*N(:,:,3)*L2m(:,:,3)*N(:,:,3)'...
            + De(4)*N(:,:,4)*L2m(:,:,4)*N(:,:,4)'...
            + De(5)*N(:,:,5)*L2m(:,:,5)*N(:,:,5)'...
            + De(6)*N(:,:,6)*L2m(:,:,6)*N(:,:,6)'; % undamaged part of stiffness operator in Voigt notation
end

%% Check stiffness/compliance operator decomposition
if ischarin('check',varargin)
    tol = 1e-12;
    D = calc_opmat(mat,elem,xnode,xgauss); % stiffness operator in Voigt notation
    decompD = max(norm(D - (Dp+Dm))/norm(D),[],'all'); if decompD>tol, decompD, end
end
