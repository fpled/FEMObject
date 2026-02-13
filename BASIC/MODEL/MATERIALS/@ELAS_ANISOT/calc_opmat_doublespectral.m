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
[Ve,De] = eig(MYDOUBLEND(D),'vector','sort');

% Eigenprojectors
% switch dim
%     case 1
%         M1 = Ve(:,1)*Ve(:,1)';
%     case 2
%         M1 = Ve(:,1)*Ve(:,1)';
%         M2 = Ve(:,2)*Ve(:,2)';
%         M3 = Ve(:,3)*Ve(:,3)';
%         if isaxi(elem)
%             M4 = Ve(:,4)*Ve(:,4)';
%         end
%     case 3
%         M1 = Ve(:,1)*Ve(:,1)';
%         M2 = Ve(:,2)*Ve(:,2)';
%         M3 = Ve(:,3)*Ve(:,3)';
%         M4 = Ve(:,4)*Ve(:,4)';
%         M5 = Ve(:,5)*Ve(:,5)';
%         M6 = Ve(:,6)*Ve(:,6)';
% end

%% Check numerical spectral (eigenvalue) decomposition of elasticity tensor
% tol = 1e-12;
% if verLessThan('matlab','9.1') % compatibility (<R2016b)
%     Dnum = Ve*diag(De)*Ve';
% else
%     Dnum = Ve*(De.*Ve');
% end
% decompDnum = max(norm(D - Dnum)/norm(D),[],'all'); if decompDnum>tol, decompDnum, end

%% Check numerical and analytical spectral (eigenvalue) decompositions of elasticity tensor
% tol = 1e-12;
% switch dim
%     case 1
%         Dnum = De(1)*M1;
%     case 2
%         if isaxi(elem)
%             Dnum = De(1)*M1 + De(2)*M2 + De(3)*M3 + De(4)*M4;
%         else
%             Dnum = De(1)*M1 + De(2)*M2 + De(3)*M3;
%         end
%     case 3
%         Dnum = De(1)*M1 + De(2)*M2 + De(3)*M3...
%             + De(4)*M4 + De(5)*M5 + De(6)*M6;
% end
% errD = max(norm(D-Dnum)/norm(D),[],'all'); if errD>tol, errD, end
% 
% switch dim
%     case 1
%         errM = max(norm(M1-eye(1)),[],'all'); if errM>tol, errM, end
%     case 2
%         if isaxi(elem)
%             errM = max(norm(M1+M2+M3+M4-eye(4)),[],'all'); if errM>tol, errM, end
%         else
%             errM = max(norm(M1+M2+M3-eye(3)),[],'all'); if errM>tol, errM, end
%         end
%     case 3
%         errM = max(norm(M1+M2+M3+M4+M5+M6-eye(6)),[],'all'); if errM>tol, errM, end
% end

%% Numerical spectral (eigenvalue) decomposition of eigenvectors
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
