function u = tt_als(A,b,u,opts)
% function u = tt_als(A,b,u,opts)
    if nargin == 3
        opts = struct()
    end
    opts = default_opts(opts);
    max_iter = opts.max_iter; % the "." is very slow in MATLAB
    stag_tol = opts.stag_tol; 
    res_tol = opts.res_tol;

    norm_b = norm(orth(b));
    eval_res = @(x) norm(orth(A*x-b))/norm_b;
    eval_stag = @(x0,x) norm(orth(x-x0))/norm(orth(x0));
    eval_energy = @(x) dot(x,0.5*A*x-b);
    
    if opts.display
        fprintf('TT_ALS\n\t\tinitial residual:\t%d\n',eval_res(u));
        fprintf('\t\tinitial energy:\t\t%d\n',eval_energy(u));
    end
    
    d = A.order;
    i = 1;
    res = 1;
    stag = 1;

    % basis is orthogonalized once
    u.space = orth(u.space); 
    % B, C and D do not change across iterations
    B = dot(A.space*u.space,u.space);
    for mu = 1:d
        sz = [u.space.dim(mu),A.space.dim(mu),u.space.dim(mu)];
        B{mu} = TUCKER_CORE(reshape(B{mu},sz),3,sz);
        % Let's pray for the right permutation
        % it works in all cases if all matrices of A.space are symmetric
        % otherwise it is B{mu} = permute(B{mu},[2 1 3]); 
        B{mu} = permute(B{mu},[2 3 1]);
    end
    C = cell(d,1);
    % the TUCKER_CORE class checks any size problems while
    % using ttt
    for mu = 1:d
        C{mu} = ttt(B{mu},A.core.a{mu},1,2);
    end
    D = dot(b.space,u.space);
    
    while (i <= max_iter) && (res > res_tol) && (stag > stag_tol)
        u0 = u;
        u = orth(u);
        uc = u.core;
        for lambda = 1:d
            uc = orth(uc,lambda); % orth w.r.t. the lambdas
                                  % dimension
            N = prod(uc.a{lambda}.rank);
            
            %% REDUCED OPERATOR
            L = TUCKER_CORE(1,3,ones(1,3));
            for mu = 1:lambda-1
                umu = uc.a{mu};
                L = ttt(L,umu,3,1); % L_{i'j'nk} = L_{i'j'k'} umu_{k'nk}
                L = ttt(L,C{mu},[1 3],[3 2]); % L_{j'kmi} = L_{i'j'nk} C_{mni'i}
                L = ttt(L,umu,[1 3],[1 2]); % L_{kij} = L_{j'kmi} umu_{j'mj}
                L = permute(L,[2 1 3]);
            end
            R = TUCKER_CORE(1,3,ones(1,3));
            for mu = d:-1:lambda+1
                umu = uc.a{mu};
                R = ttt(R,umu,2,3); % R_{i'k'jm} = R_{i'j'k'} umu_{jmj'}
                R = ttt(R,C{mu},[1 4],[4 1]); % R_{k'jni} = R_{i'k'jm} C_{mnii'}
                R = ttt(R,umu,[1 3],[3 2]); % R_{jik} = R_{k'jni} umu_{knk'}
                R = permute(R,[2 1 3]);
            end
            Ared = ttt(L,C{lambda},1,3); % Ared_{imjnt} = L_{sim} C_{jnst}
            Ared = ttt(Ared,R,5,1); % Ared_{imjnko} = Ared_{imjnt} R_{tko}
            Ared = permute(Ared,[1 3 5 2 4 6]);
            Ared = reshape(Ared.a,N,N);

            %% REDUCED RHS
            E = 1;
            for mu = 1:lambda-1
                E = ttm(b.core.a{mu},E',1); % E_{j'ki} = b_{i'ki} E_{i'j'}
                EE = ttm(uc.a{mu},D{mu},2); % EE_{j'kj} = umu_{j'lj} D_{kl}
                E = ttt(E,EE,[1 2],[1 2]); % E_{ij} = E_{j'ki} EE_{j'kj}
                E = E.a;
            end
            F = 1;
            for mu = d:-1:lambda+1
                F = ttm(b.core.a{mu},F',3); % F_{ikj'} = b_{iki'} F_{i'j'}
                F = ttm(F,D{mu}',2); % F_{ilj'} = F_{ikj'} D_{kl}
                F = ttt(F,uc.a{mu},[2 3],[2 3]); % F_{ij} = F_{ilj'} umu_{jlj'}
                F = F.a;
            end
            bred = b.core.a{lambda};
            bred = ttm(bred,{E',D{lambda}',F'},1:3);
            bred = reshape(bred.a,N,1);
            
            ulambda = Ared\bred;
            
            uc.a{lambda}.a(:) = ulambda(:);
        end
        u = LRTENSOR(uc,u.space);
        res = eval_res(u);
        stag = eval_stag(u0,u);
        energy = dot(u,0.5*A*u-b);
        if opts.display
            fprintf('iter %d -- stag %d -- res %d -- energy %d\n',...
                    i,stag,res,energy);
        end
        i = i + 1;
    end
end

function opts = default_opts(opts)
    if ~isfield(opts,'max_iter'); opts.max_iter = 10; end
    if ~isfield(opts,'stag_tol'); opts.stag_tol = 1e-6; end
    if ~isfield(opts,'res_tol'); opts.res_tol = 1e-6; end
    if ~isfield(opts,'display'); opts.display = true; end
end

