function mov = evol_sol_cell(T,f,varargin)
% function mov = evol_sol_cell(T,f,varargin)
% T : TIMEMODEL
% f : cell of functions
% varargin : arguments for function plot_sol
% 'rescale' : fixed scale (true) or automatically adapted scale (false by default)
% mov = movie

rescale = getcharin('rescale',varargin,false);
interfaces = getcharin('interfaces',varargin,[]);

[t,rep]=gettevol(T);

if rescale
    
    qtmax = cell(1,length(f));
    qtmin = cell(1,length(f));
    for i=1:length(t)
        if ischarin('sigma',varargin)
            q = calc_sigma(varargin{1}{i},f{rep(i)},varargin{2:end});
            k = getcharin('sigma',varargin);
            q = sigmacompo(q,k,varargin{1}{i});
        elseif ischarin('epsilon',varargin)
            q = calc_epsilon(varargin{1}{i},f{rep(i)},varargin{2:end});
            k = getcharin('epsilon',varargin);
            q = sigmacompo(q,k,varargin{1}{i});
        elseif ischarin('energyint',varargin)
            q = calc_energyint(varargin{1}{i},f{rep(i)},varargin{2:end});
        elseif ischarin('displ',varargin)
            q = unfreevector(varargin{1}{i},f{rep(i)});
            k = getcharin('displ',varargin);
            switch k
                case 1
                    numddl = findddl(varargin{1}{i},'UX');
                case 2
                    numddl = findddl(varargin{1}{i},'UY');
                case 3
                    numddl = findddl(varargin{1}{i},'UZ');
            end
            q = q(numddl);
        elseif ischarin('rotation',varargin)
            q = unfreevector(varargin{1}{i},f{rep(i)});
            k = getcharin('rotation',varargin);
            switch k
                case 1
                    numddl = findddl(varargin{1}{i},'RX');
                case 2
                    numddl = findddl(varargin{1}{i},'RY');
                case 3
                    numddl = findddl(varargin{1}{i},'RZ');
            end
            q = q(numddl);
        else
            q = f{rep(i)};
        end
        
        qtmax{i} = max(q);
        qtmin{i} = min(q);
    end
    qtmax = max([qtmax{:}]);
    qtmin = min([qtmin{:}]);
    
    T = setevolparam(T,'caxis',[qtmin,qtmax]);
    if strfind(rescale,'z')
        T = setevolparam(T,'zlim',[qtmin,qtmax]);
    elseif strfind(rescale,'y')
        T = setevolparam(T,'ylim',[qtmin,qtmax]);
    end
    varargin = delcharin('rescale',varargin);
    
end

fontsize = getevolparam(T,'FontSize');
plotiter = getevolparam(T,'plotiter');
plottime = getevolparam(T,'plottime');
plotstep = getevolparam(T,'plotstep');
compt = 0;
for i=1:plotstep:length(t)
    compt=compt+1;
    execute(T.evolparam,'before');
    
    plot_sol(f{rep(i)},varargin{1}{i},varargin{2:end});
    
    if plotiter
        title(['iter ' num2str(i,'%d')],'FontSize',fontsize)
    elseif plottime
        title(['time ' num2str(t(i),'%.2f') ' s'],'FontSize',fontsize)
    end
    
    execute(T.evolparam,'after')
    
    if nargout>=1
        mov(compt) = getframe(gcf);
    end
end

function se = sigmacompo(se,ksigma,S)

if isa(ksigma,'char') && strcmp(ksigma,'mises')
    switch getindim(S)
        case 1
            se = se(1);
        case 2
            se1 = se(1);
            se2 = se(2);
            se3 = se(3);
            tracese = 1/3*(se1+se2);
            se1 = se1 - tracese;
            se2 = se2 - tracese;
            se = sqrt(3/2*(se1.*se1 + se2.*se2 + 2*se3.*se3));
        case 3
            se1 = se(1);
            se2 = se(2);
            se3 = se(3);
            se4 = se(4);
            se5 = se(5);
            se6 = se(6);
            tracese = 1/3*(se1+se2+se3);
            se1 = se1 - tracese;
            se2 = se2 - tracese;
            se3 = se3 - tracese;
            se = sqrt(3/2*(se1.*se1 + se2.*se2 + se3.*se3...
                + 2*(se4.*se4 + se5.*se5 + se6.*se6)));
    end
else
    se = se(ksigma);
end

return
