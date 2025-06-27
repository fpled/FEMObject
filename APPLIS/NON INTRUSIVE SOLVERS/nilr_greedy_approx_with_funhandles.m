function [V,Lambda,result] = nilr_greedy_approx_with_funhandles(pw_J,pw_R,sz_v,opts)
% DO NOT USE (NICE BUT SLOW)
% Functional implementation of the greedy algorithm
%

opts = defaultopts(opts);

result.time = zeros(1,opts.maxcorrec);
result.error = result.time;
result.results = cell(1,opts.maxcorrec);

[V,Lambda,result.results{1}] = nilr_altern_mini(pw_J,pw_R,sz_v,opts);
result.time(1) = result.results{1}.time;
ucur = PCMATRIX(V*Lambda',[sz_v 1],opts.PC);

for k = 2:opts.maxcorrec
    u0 = ucur;
    pw_Jk = @(p,du) pw_J(p,legendre_approx(ucur,ind,p)+du);
    pw_Rk = @(p,du) pw_R(p,legendre_approx(ucur,ind,p)+du);
    [dV,dLambda,result.results{k}] = nilr_altern_mini(pw_Jk,pw_Rk,sz_v,opts);
    V = [V dV];
    Lambda = [Lambda dLambda];

    ucur = PCMATRIX(V*Lambda',[sz_v 1],opts.PC);
    
    result.error(k) = norm(ucur-u0)/norm(u0);
    result.time(k) = result.results{k}.time;
end

end

function opts = defaultopts(opts)
    if ~isfield(opts,'maxcorrec');opts.maxcorrec = 10;end
    if ~isfield(opts,'PC');error('Need a POLYCHAOS');end
end