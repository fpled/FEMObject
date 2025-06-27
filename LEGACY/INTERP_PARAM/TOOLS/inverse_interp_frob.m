function [P,residual,out] = inverse_interp_frob(A,Y,V,opts)
%  function P = inverse_interp_frob(A,Y,V)
% 
% Frobenius norm projection : min_{P\in Y} || I-PA ||_F
% 
% If V is specified : min_{P\in Y} || (I-PA)V ||_F

if nargin <=2 || isempty(V)
    V=eye(A.sz(1,:));
end
if nargin <=3 || isempty(opts)
    opts = struct();
end
opts = default_opts(opts);

%% M_{ij}(\xi)

opts_eim.max_basis=1000;
opts_eim.tol=1e-10;

PhiA = mtimes(A.space,A.space,2);
PhiA = cellfun(@spdiags,PhiA.u{1},'UniformOutput',0);
PhiA = [PhiA{:}];
[~,~,PivotElementM_1,PivotElementM_2]=eim(PhiA',opts_eim);
M_Phi = PhiA(:,PivotElementM_2) *inv( PhiA(PivotElementM_1,PivotElementM_2) );
clear PhiA

n_evalM = length(PivotElementM_1);


if nargout == 3
    nargout3=1;
    outPAV_M=cell(n_evalM,1);
    out.PivotElementM_1=PivotElementM_1;
    out.M_Phi=M_Phi;
else
    nargout3=0
end
M_xi=cell(n_evalM,1);

parfor i=1:n_evalM
    if opts.display
        fprintf('\n%d/%d ', i,n_evalM)
    end
    
    AV=eval_param_operator(A,PivotElementM_1(i))*V;
    PAV = cellfun(@(y) y*AV,Y, 'UniformOutput',0);
    PAV = cellfun(@(x) x(:),PAV, 'UniformOutput',0);
    PAV = [PAV{:}];
    if nargout3
        outPAV_M{i}=PAV;
    end
    M_xi{i} = PAV'*PAV;
    
    
end

if nargout == 3
    out.PAV_M=outPAV_M;
    out.M_xi=M_xi;
end

%% S_i(\xi)


PhiA=cellfun(@spdiags,A.space.u{2},'UniformOutput',0);
PhiA=[PhiA{:}];
[~,~,PivotElementS_1,PivotElementS_2]=eim(PhiA',opts_eim);
S_Phi = PhiA(:,PivotElementS_2) / PhiA(PivotElementS_1,PivotElementS_2);
clear PhiA

n_evalS = length(PivotElementS_1);

if nargout == 3
    outPAV_S=cell(n_evalS,1);
    out.PivotElementS_1=PivotElementS_1;
    out.S_Phi=S_Phi;
end
S_xi=cell(n_evalS,1);
parfor i=1:n_evalS
    if opts.display
        fprintf('\n%d/%d ', i,n_evalS)
    end
    
    AV=eval_param_operator(A,PivotElementS_1(i))*V;
    S_xi{i} = cellfun(@(y) trace( V'*(y*AV) ),Y);
%     if nargout == 3
%         outPAV_M{i}=PAV;
%     end
    
    
end

if nargout == 3
    out.PAV_M=outPAV_M;
    out.S_xi=S_xi;
end

%%

[Phi,residual]=MS_solve(M_xi,M_Phi,S_xi,S_Phi);

PT = TSPACE_OPERATORS( [{Y};{param_operator( Phi )}] );
PC = CANONICAL_CORE( ones(1,PT.dim(1)),2 );
P = LRTENSOR(PC,PT);


residual=residual + trace(V'*V);
residual(residual<0)=0;



end

function opts = default_opts(opts)
    if ~isfield(opts,'display'); opts.display= true; end
    if ~isfield(opts,'tol'); opts.tol = 1e-12; end
    
%     if ~isfield(opts,'orth'); opts.orth = true; end
%     if ~isfield(opts,'metric'); opts.metric = @(x) norm(x,Inf); end
%     if ~isfield(opts,'dot'); opts.dot = []; end
end




