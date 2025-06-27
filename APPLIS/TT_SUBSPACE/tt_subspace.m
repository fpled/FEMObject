function u = tt_subspace(A,b,opts)
% function u = tt_subspace(A,b,opts)
    
    if nargin == 2
        opts = struct();
    end
    opts = defaultopts(opts);
    
    optim_core = opts.optim_core;
    opts_optim_core = opts.opts_optim_core;
    opts_rankone = opts.opts_rankone;
    
    max_iter = opts.max_iter;
    stag_tol = opts.stag_tol;
    res_tol = opts.res_tol;
    
    opts_rankone.display = false;
    opts_optim_core.display = false;
    
    res = 1;
    stag = 1;
    energy = 1;
    
    norm_b = norm(orth(b));
    
    eval_stag = @(x0,x) norm(orth(x-x0))/norm(orth(x));
    eval_res = @(x) norm(orth(A*x-b))/norm_b;
    eval_energy = @(x) dot(x,0.5*A*x-b);
    
    u = rankoneapproxsolALS(A,b,opts_rankone);
    u = LRTENSOR(TT_CORE(u.core),u.space); % conversion to TT_CORE
    
    res = eval_res(u);
    energy = eval_energy(u);
    if opts.display
        fprintf('Iter 1 -- stag %d -- res %d -- energy %d\n', ...
                stag, res, energy);
    end
    i = 2;
    while (i <= max_iter) && (stag >= stag_tol) && (res >= res_tol)
        u0 = u;
        ru = orth(b-A*u);
        v = rankoneapproxsolALS(A,ru,opts_rankone);
        % dirty subspace enrichment
        u = orth(u + v);
        
        % computation of the coordinates
        if strcmp(optim_core,'als')
            u = tt_als(A,b,u,opts_optim_core);
        elseif strcmp(optim_core,'dmrg')
            u = tt_dmrg(A,b,u,opts_optim_core);
        end
        
        stag = eval_stag(u0,u);
        res = eval_res(u);
        energy = eval_energy(u);
        if opts.display
            fprintf('Iter %d -- stag %d -- res %d -- energy %d\n', ...
                    i, stag, res, energy);
        end
        i = i+1;
    end
end

function opts = defaultopts(opts)
    if ~isfield(opts,'optim_core'); opts.optim_core = 'als'; end
    if ~isfield(opts,'opts_optim_core'); opts.opts_optim_core = struct(); end
    if ~isfield(opts,'opts_rankone'); opts.opts_rankone = struct(); end
    if ~isfield(opts,'max_iter'); opts.max_iter = 10; end
    if ~isfield(opts,'stag_tol'); opts.stag_tol = 1e-5; end
    if ~isfield(opts,'res_tol'); opts.res_tol = 1e-5; end
    if ~isfield(opts,'display'); opts.display = true; end
end
