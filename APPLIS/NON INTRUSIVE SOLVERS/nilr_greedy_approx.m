function [V,Lambda,result] = nilr_greedy_approx(pw_J,pw_R,sz_v,opts)
% function [V,Lambda,result] = nilr_greedy_approx(pw_J,pw_R,sz_v,opts)
%
% Compute a low-rank approximation V*Lambda'of the solution u in a greedy
% fashion satisfying 
% pw_J(p,u(p)) = min pw_J(p,v), for all p.
% 
% pw_R is the residual of the equation, that is 
% pw_R(p,v) = -\nabla_v pw_J(p,v)
% 
% **** Option concerning the greedy approximation is
% maxcorrec : number of successive corrections applied to the approximation
%
% For other options see nilr_altern_mini / nilr_altern_mini_for_correc
%

opts = defaultopts(opts);

result.time = zeros(1,opts.maxcorrec);
result.error = result.time;
result.results = cell(1,opts.maxcorrec);

[V,Lambda,result.results{1}] = nilr_altern_mini(pw_J,pw_R,sz_v,opts);
result.time(1) = result.results{1}.time;
ucur = V*Lambda';

for k = 2:opts.maxcorrec
    u0 = ucur;
    
    [dV,dLambda,result.results{k}] = nilr_altern_mini_for_correc(pw_J,pw_R,sz_v,...
        V,Lambda,opts);
    V = [V dV];
    Lambda = [Lambda dLambda];
    ucur = V*Lambda';
    
    result.error(k) = norm(ucur-u0,'fro')/norm(u0,'fro');
    result.time(k) = result.results{k}.time;
end

end

function opts = defaultopts(opts)
    if ~isfield(opts,'maxcorrec');opts.maxcorrec = 10;end
end