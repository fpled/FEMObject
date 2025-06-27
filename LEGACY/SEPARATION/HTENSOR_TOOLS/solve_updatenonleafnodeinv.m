function [P] = solve_updatenonleafnodeinv(A,P,t,opts)
% function [P] = solve_updatenonleafnodeinv(A,P,t,opts)

% N=numel(u.B{t});
% uAu=zeros(N);
% Au=apply_mat_to_vec(A,u);
%
% sz_A = size(A.B{t});
% sz_A(end+1:3) = 1;
% sz_u = size(u.B{t});
% sz_u(end+1:3) = 1;
%
% M = reduceinnerprodatnode(u,Au,t);
% Mb = reduceinnerprodatnode(u,b,t);

if nargin == 3
    opts.residual = 0;
else
    if ~isfield(opts,'residual');opts.residual=0;end
end

N = numel(P.B{t});
PCP = zeros(N);

n = sqrt(size(A));
I = sepmatrixtohtensor(sepopereye(n));

if opts.residual
    B = transpose_htensor(A);
    C = apply_mat_to_mat(A,B,n);
else
    B = I;
    C = A;
end

CP = apply_mat_to_mat(C,P,n);
sz_C = size(C.B{t});
sz_C(end+1:3) = 1;
sz_P = size(P.B{t});
sz_P(end+1:3) = 1;

M = reduceinnerprodatnode(P,CP,t);
BM = reduceinnerprodatnode(P,B,t);

if t==1 %t is root node
    %     for n=1:N
    %         Br=zeros(size(u.B{1}));
    %         Br(n)=1;
    %         AuB=Au.B{1};
    %         AuB(:, :) = kron(A.B{1}(:, :), Br(:, :));
    %         B = M{1}*AuB*M{2}';
    %         uAu(:,n)=B(:);
    %     end
    %     ub = ttm(b.B{1}, Mb(1:2), 1:2);
    %     ub=ub(:);
    for n=1:N
        Br=zeros(size(P.B{1}));
        Br(n)=1;
        CPB=CP.B{1};
        CPB(:, :) = kron(C.B{1}(:, :), Br(:, :));
        MBM = M{1}*CPB*M{2}';
        PCP(:,n)=MBM(:);
    end
    PB = ttm(B.B{1}, BM(1:2), 1:2);
    PB=PB(:);
else
    %     Bb=ttm(b.B{t},Mb,[1 2 3]);
    %     szuBt=size(u.B{t}); %Suppressing "." make things way faster
    %     szAuB=size(Au.B{t});
    %     AuB=zeros(szAuB);
    %     iA=((1:sz_A(1)) - 1) * sz_u(1);
    %     jA=((1:sz_A(2)) - 1) * sz_u(2);
    %     kA=((1:sz_A(3)) - 1) * sz_u(3);
    %     ABt=A.B{t};
    %     for n=1:N
    %         if n>1
    %             % Allocating at each iteration takes too much time,
    %             % almost as much as "B = ttm(AuB,M,[1 2 3])".
    %             % Erasing previous entries is much faster.
    %             AuB(iAuB,jAuB,kAuB)=zeros(numel(iAuB),numel(jAuB),numel(kAuB));
    %         end
    %         [iu,ju,ku] = ind2sub(szuBt,n);
    %         iAuB=iu+iA;
    %         jAuB=ju+jA;
    %         kAuB=ku+kA;
    %         AuB(iAuB,jAuB,kAuB)=ABt;
    %         B = ttm(AuB,M,[1 2 3]);
    %         uAu(n,:) = B(:);
    %     end
    %     ub = Bb(:);
    BB = ttm(B.B{t},BM,[1 2 3]);
    szPBt=size(P.B{t}); %Suppressing "." make things way faster
    szCPB=size(CP.B{t});
    CPB=zeros(szCPB);
    iC=((1:sz_C(1)) - 1) * sz_P(1);
    jC=((1:sz_C(2)) - 1) * sz_P(2);
    kC=((1:sz_C(3)) - 1) * sz_P(3);
    CBt=C.B{t};
    for n=1:N
        if n>1
            % Allocating at each iteration takes too much time,
            % almost as much as "B = ttm(AuB,M,[1 2 3])".
            % Erasing previous entries is much faster.
            CPB(iCPB,jCPB,kCPB)=zeros(numel(iCPB),numel(jCPB),numel(kCPB));
        end
        [iP,jP,kP] = ind2sub(szPBt,n);
        iCPB=iP+iC;
        jCPB=jP+jC;
        kCPB=kP+kC;
        CPB(iCPB,jCPB,kCPB)=CBt;
        MBM = ttm(CPB,M,[1 2 3]);
        PCP(n,:) = MBM(:);
    end
    PB = BB(:);
end

% u.B{t} = reshape(uAu\ub,size(u.B{t}));
P.B{t} = reshape(PCP\PB,size(P.B{t}));

end

