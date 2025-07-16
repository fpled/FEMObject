function mov = evol_sol(T,f,varargin)
% function mov = evol_sol(T,f,varargin)
% T : TIMEMODEL
% f : function or cell of functions
% varargin : arguments for function plot_sol
% 'rescale' : fixed scale (true) or automatically adapted scale (false by default)
% mov = movie

rescale = getcharin('rescale',varargin,false);
interfaces = getcharin('interfaces',varargin,[]);

[t,rep]=gettevol(T);
if ~iscell(f)
    f = {f};
end
for j=1:length(f)
    if ~isa(f{j},'TIMEMATRIX') && ~isa(f{j},'TIMERADIALMATRIX')
        f{j} = TIMEMATRIX(f{j},T);
    end
end
if ~iscell(varargin{1})
    varargin{1} = varargin(1);
end

if rescale
    
    qtmax = cell(1,length(f));
    qtmin = cell(1,length(f));
    for j=1:length(f)
        if ischarin('sigma',varargin)
            q = calc_sigma(varargin{1}{j},f{j},varargin{2:end});
            k = getcharin('sigma',varargin);
            q = sigmacompo(q,k,varargin{1}{j});
        elseif ischarin('epsilon',varargin)
            q = calc_epsilon(varargin{1}{j},f{j},varargin{2:end});
            k = getcharin('epsilon',varargin);
            q = sigmacompo(q,k,varargin{1}{j});
        elseif ischarin('energyint',varargin)
            q = calc_energyint(varargin{1}{j},f{j},varargin{2:end});
        elseif ischarin('displ',varargin)
            q = unfreevector(varargin{1}{j},f{j});
            k = getcharin('displ',varargin);
            switch k
                case 1
                    numddl = findddl(varargin{1}{j},'UX');
                case 2
                    numddl = findddl(varargin{1}{j},'UY');
                case 3
                    numddl = findddl(varargin{1}{j},'UZ');
            end
            q = getcompo(q,numddl);
        elseif ischarin('rotation',varargin)
            q = unfreevector(varargin{1}{j},f{j});
            k = getcharin('rotation',varargin);
            switch k
                case 1
                    numddl = findddl(varargin{1}{j},'RX');
                case 2
                    numddl = findddl(varargin{1}{j},'RY');
                case 3
                    numddl = findddl(varargin{1}{j},'RZ');
            end
            q = getcompo(q,numddl);
        else
            q = f{j};
        end
        
        if isa(q,'TIMEMATRIX') && ~israndom(q)
            qtmax{j} = max(max(q));
            qtmin{j} = min(min(q));
        elseif isa(q,'TIMERADIALMATRIX')
            qtmax{j}=-Inf;
            qtmin{j}=Inf;
            for kk=1:max(1,floor(length(rep)/5)):length(rep)
                qtmax{j} = max(qtmax{j},max(max(getmatrixatstep(q,rep(kk)))));
                qtmin{j} = min(qtmin{j},min(min(getmatrixatstep(q,rep(kk)))));
            end
        end
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
    
    for j=1:length(f)
        if ~ismember(j,interfaces)
            plot_sol(getmatrixatstep(f{j},rep(i)),varargin{1}{j},varargin{2:end});
        elseif ~ischarin('sigma',varargin) && ~ischarin('epsilon',varargin) && ~ischarin('energyint',varargin)
            plot_sol(getmatrixatstep(f{j},rep(i)),varargin{1}{j},'FaceColor','none','EdgeColor','k',varargin{2:end});
        end
        hold on
    end
    hold off
    
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
            se = getcompo(se,1);
        case 2
            se1 = getcompo(se,1);
            se2 = getcompo(se,2);
            se3 = getcompo(se,3);
            tracese = 1/3*(se1+se2);
            se1 = se1 - tracese;
            se2 = se2 - tracese;
            T = gettimemodel(se);
            se = sqrt(3/2*(getvalue(se1).*getvalue(se1) + getvalue(se2).*getvalue(se2) + 2*getvalue(se3).*getvalue(se3)));
            se = TIMEMATRIX(se,T);
        case 3
            se1 = getcompo(se,1);
            se2 = getcompo(se,2);
            se3 = getcompo(se,3);
            se4 = getcompo(se,4);
            se5 = getcompo(se,5);
            se6 = getcompo(se,6);
            tracese = 1/3*(se1+se2+se3);
            se1 = se1 - tracese;
            se2 = se2 - tracese;
            se3 = se3 - tracese;
            T = gettimemodel(se);
            se = sqrt(3/2*(getvalue(se1).*getvalue(se1) + getvalue(se2).*getvalue(se2) + getvalue(se3).*getvalue(se3)...
                + 2*(getvalue(se4).*getvalue(se4) + getvalue(se5).*getvalue(se5) + getvalue(se6).*getvalue(se6))));
            se = TIMEMATRIX(se,T);
    end
else
    se = getcompo(se,ksigma);
end

return
