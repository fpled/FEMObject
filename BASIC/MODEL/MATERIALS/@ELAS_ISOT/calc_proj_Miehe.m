function [Pp,Pm] = calc_proj_Miehe(mat,elem,xnode,xgauss,se,varargin)
% function [Pp,Pm] = calc_proj_Miehe(mat,elem,xnode,xgauss,se)

dim = getdim(elem);
split = getparam(mat,'PFS'); % phase field split

if strcmpi(split,'stress')
    P = calc_proj_notation(elem);
    B = P*P';
    se = B*se; % stress tensor converted as strain tensor in Voigt notation
end

%% Strain/Stress tensor
switch dim
    case 1
        Se = se(1);
    case 2
        Se = [se(1) se(3)/2;
            se(3)/2 se(2)];
    case 3
        Se = [se(1) se(4)/2 se(5)/2;
            se(4)/2 se(2) se(6)/2;
            se(5)/2 se(6)/2 se(3)];
end

%% Numerical spectral (eigenvalue) decomposition of strain/stress tensor
% % Eigenvalues and eigenvectors
% [vecnum,valnum] = eig(Se,'vector','sort');
% 
% % Eigenprojectors
% switch dim
%     case 1
%         M1num = vecnum(:,1)*vecnum(:,1)';
%     case 2
%         M1num = vecnum(:,1)*vecnum(:,1)';
%         M2num = vecnum(:,2)*vecnum(:,2)';
%     case 3
%         M1num = vecnum(:,1)*vecnum(:,1)';
%         M2num = vecnum(:,2)*vecnum(:,2)';
%         M3num = vecnum(:,3)*vecnum(:,3)';
% end

%% Check numerical spectral (eigenvalue) decomposition of strain/stress tensor
% tol = 1e-12;
% if verLessThan('matlab','9.1') % compatibility (<R2016b)
%     Senum = vecnum*diag(valnum)*vecnum';
% else
%     Senum = vecnum*(valnum.*vecnum');
% end
% decompSenum = max(norm(Se - Senum)/norm(Se),[],'all'); if decompSenum>tol, decompSenum, end

%% Analytical spectral (eigenvalue) decomposition of strain/stress tensor
switch dim
    case 1
        M1 = zerosND(size(Se));
        M1(1,1) = 1;
    case 2
        % Invariants
        I1 = trace(Se);
        I2 = det(Se);
        
        % Eigenvalues (in ascending order)
        val = [(I1 - sqrt(I1^2 - 4*I2))/2, (I1 + sqrt(I1^2 - 4*I2))/2]';
        k = find(ne(val(1),val(2))); % positions of 2 distinct eigenvalues
        
        % Eigenprojectors
        M1 = zerosND(size(Se));
        M1(1,1) = 1;
        M1(:,:,k) = (Se(:,:,k) - val(2,:,k)*eye(2))/(val(1,:,k)-val(2,:,k));
        % M2 = zerosND(size(Se));
        % M2(2,2) = 1;
        % M2(:,:,k) = (Se(:,:,k) - val(1,:,k)*eye(2))/(val(2,:,k)-val(1,:,k));
        M2 = eye(2) - M1;
        
    case 3
        % Invariants
        I1 = trace(Se);
        I2 = (trace(Se)^2 - trace(Se*Se))/2;
        I3 = det(Se);
        
        % Eigenvalues (in ascending order)
        val = I1/3*ones(3,1); % initialization : case of 3 equal eigenvalues
        g = I1^2 - 3*I2;
        arg = (2*I1^3 -9*I1*I2 +27*I3)/2/g^(3/2); % -1 <= arg <= 1
        theta = acos(arg)/3; % Lode's angle such that 0 <= theta <= pi/3
        kmin = find(arg==1); % positions of double minimum eigenvalue
        kmax = find(arg==-1); % positions of double maximum eigenvalue
        k = find(ne(g,0));
        k = setdiff(k,kmin);
        k = setdiff(k,kmax); % positions of 3 distinct eigenvalues
        % Case of double minimum eigenvalue
        val(:,:,kmin) = val(:,:,kmin) + sqrt(g(:,:,kmin))*[-1 -1 2]'/3;
        % Case of double maximum eigenvalue
        val(:,:,kmax) = val(:,:,kmax) + sqrt(g(:,:,kmax))*[-2 1 1]'/3;
        % Case of 3 distinct eigenvalues
        val(:,:,k) = val(:,:,k) + sqrt(g(:,:,k))*...
            [cos(2*pi/3+theta(:,:,k)); cos(2*pi/3-theta(:,:,k)); cos(theta(:,:,k))]*2/3;
        
        % Eigenprojectors
        M1 = zerosND(size(Se));
        M1(1,1) = 1;
        M3 = zerosND(size(Se));
        M3(3,3) = 1;
        % Case of double minimum eigenvalue
        M3(:,:,kmin) = (Se(:,:,kmin) - (I1(:,:,kmin)-sqrt(g(:,:,kmin)))*eye(3)/3)./sqrt(g(:,:,kmin));
        M1(:,:,kmin) = (eye(3) - M3(:,:,kmin))/2;
        % Case of double maximum eigenvalue
        M1(:,:,kmax) = ((I1(:,:,kmax)-sqrt(g(:,:,kmax)))*eye(3)/3 - Se(:,:,kmax))./sqrt(g(:,:,kmax));
        M3(:,:,kmax) = (eye(3) - M1(:,:,kmax))/2;
        % Case of 3 distinct eigenvalues
        M1(:,:,k) = (Se(:,:,k)-val(2,:,k)*eye(3))/(val(1,:,k)-val(2,:,k))*...
            (Se(:,:,k)-val(3,:,k)*eye(3))/(val(1,:,k)-val(3,:,k));
        M3(:,:,k) = (Se(:,:,k)-val(1,:,k)*eye(3))/(val(3,:,k)-val(1,:,k))*...
            (Se(:,:,k)-val(2,:,k)*eye(3))/(val(3,:,k)-val(2,:,k));
        % M2 = zerosND(size(Se));
        % M2(2,2) = 1;
        % M2(:,:,kmin) = M1(:,:,kmin);
        % M2(:,:,kmax) = M3(:,:,kmax);
        % M2(:,:,k) = (Se(:,:,k)-val(1,:,k)*eye(3))/(val(2,:,k)-val(1,:,k))*...
        %     (Se(:,:,k)-val(3,:,k)*eye(3))/(val(2,:,k)-val(3,:,k));
        M2 = eye(3) - (M1 + M3);
        
end

%% Check numerical and analytical spectral (eigenvalue) decompositions of strain/stress tensor
% tol = 1e-12;
% switch dim
%     case 1
%         Senum = valnum(1)*M1num;
%         Seana = val(1)*M1;
%     case 2
%         Senum = valnum(1)*M1num + valnum(2)*M2num;
%         Seana = val(1)*M1 + val(2)*M2;
%     case 3
%         Senum = valnum(1)*M1num + valnum(2)*M2num + valnum(3)*M3num;
%         Seana = val(1)*M1 + val(2)*M2 + val(3)*M3;
% end
% errSenum = max(norm(Se-Senum)/norm(Se),[],'all'); if errSenum>tol, errSenum, end
% errSeana = max(norm(Se-Seana)/norm(Se),[],'all'); if errSeana>tol, errSeana, end
% 
% k = find(abs(val)>tol);
% normval = double(abs(val));
% errval = double(abs(val-valnum));
% errval(k) = errval(k)./normval(k); errval = max(errval,[],'all'); if errval>tol, errval, end
% errM1 = max(norm(M1-M1num)/norm(M1),[],'all'); if errM1>tol, errM1, end
% if dim>=2
%     errM2 = max(norm(M2-M2num)/norm(M2),[],'all'); if errM2>tol, errM2, end
% end
% if dim==3
%     errM3 = max(norm(M3-M3num)/norm(M3),[],'all'); if errM3>tol, errM3, end
% end
% switch dim
%     case 1
%         errM = max(norm(M1-eye(1)),[],'all'); if errM>tol, errM, end
%         errMnum = max(norm(M1num-eye(1)),[],'all'); if errMnum>tol, errMnum, end
%     case 2
%         errM = max(norm(M1+M2-eye(2)),[],'all'); if errM>tol, errM, end
%         errMnum = max(norm(M1num+M2num-eye(2)),[],'all'); if errMnum>tol, errMnum, end
%     case 3
%         errM = max(norm(M1+M2+M3-eye(3)),[],'all'); if errM>tol, errM, end
%         errMnum = max(norm(M1num+M2num+M3num-eye(3)),[],'all'); if errMnum>tol, errMnum, end
% end

%% Computation of projectors of strain/stress tensor onto positive (extensive) and negative (compressive) parts
valp = (val+abs(val))./2;
valm = (val-abs(val))./2;

dvalp = heaviside(val);
dvalm = heaviside(-val);

% Projectors
Pp = zerosND([size(se,1),size(se,1),sizeND(se)]);
Pm = zerosND([size(se,1),size(se,1),sizeND(se)]);
switch dim
    case 2
        m1 = [M1(1,1) M1(2,2) M1(1,2)]';
        m2 = [M2(1,1) M2(2,2) M2(1,2)]';
        k = find(ne(val(1),val(2)));
        betap = dvalp(1);
        betam = dvalm(1);
        betap(:,:,k) = (valp(1,:,k)-valp(2,:,k))./(val(1,:,k)-val(2,:,k));
        betam(:,:,k) = (valm(1,:,k)-valm(2,:,k))./(val(1,:,k)-val(2,:,k));
        gammap = dvalp-betap;
        gammam = dvalm-betam;
        
        % P = calc_proj_notation(elem);
        % Id = P'\eye(3)/P;
        Id = diag([1,1,1/2]);
        M11 = m1*m1';
        M22 = m2*m2';
        Pp = betap*Id + gammap(1)*M11 + gammap(2)*M22;
        Pm = betam*Id + gammam(1)*M11 + gammam(2)*M22;
        
    case 3
        m1 = [M1(1,1) M1(2,2) M1(3,3) M1(1,2) M1(1,3) M1(2,3)]';
        m2 = [M2(1,1) M2(2,2) M2(3,3) M2(1,2) M2(1,3) M2(2,3)]';
        m3 = [M3(1,1) M3(2,2) M3(3,3) M3(1,2) M3(1,3) M3(2,3)]';
        thetap = dvalp/2;
        thetam = dvalm/2;
        k = find(ne(val(1),val(2)));
        thetap(1,:,k) = (valp(1,:,k)-valp(2,:,k))/(val(1,:,k)-val(2,:,k))/2;
        thetam(1,:,k) = (valm(1,:,k)-valm(2,:,k))/(val(1,:,k)-val(2,:,k))/2;
        k = find(ne(val(1),val(3)));
        thetap(2,:,k) = (valp(1,:,k)-valp(3,:,k))/(val(1,:,k)-val(3,:,k))/2;
        thetam(2,:,k) = (valm(1,:,k)-valm(3,:,k))/(val(1,:,k)-val(3,:,k))/2;
        k = find(ne(val(2),val(3)));
        thetap(3,:,k) = (valp(2,:,k)-valp(3,:,k))/(val(2,:,k)-val(3,:,k))/2;
        thetam(3,:,k) = (valm(2,:,k)-valm(3,:,k))/(val(2,:,k)-val(3,:,k))/2;
        
        M11 = m1*m1';
        M22 = m2*m2';
        M33 = m3*m3';
        G12 = construct_G(m1,m2);
        G13 = construct_G(m1,m3);
        G23 = construct_G(m2,m3);
        
        Pp = dvalp(1)*M11 + dvalp(2)*M22 + dvalp(3)*M33 + thetap(1)*G12 + thetap(2)*G13 + thetap(3)*G23;
        Pm = dvalm(1)*M11 + dvalm(2)*M22 + dvalm(3)*M33 + thetam(1)*G12 + thetam(2)*G13 + thetam(3)*G23;
        
end

%% Check orthogonality condition and strain/stress decomposition
% tol = 1e-12;
% P = calc_proj_notation(elem);
% sep = P*(Pp*se); % positive part of strain/stress tensor in Kelvin-Mandel notation
% sem = P*(Pm*se); % negative part of strain/stress tensor in Kelvin-Mandel notation
% se = P\se; % strain/stress tensor in Kelvin-Mandel notation
% orthpm = max(abs(sep'*sem)/abs(se'*se),[],'all'); if orthpm>tol, orthpm, end
% orthmp = max(abs(sem'*sep)/abs(se'*se),[],'all'); if orthmp>tol, orthmp, end
% decomp = max(norm(se - (sep+sem))/norm(se),[],'all'); if decomp>tol, decomp, end

end
