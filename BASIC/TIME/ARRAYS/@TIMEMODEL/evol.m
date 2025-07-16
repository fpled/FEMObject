function mov = evol(T,f,varargin)
% function mov = evol(T,f,varargin)
% T : TIMEMODEL
% f : function or cell of functions
% varargin : arguments for function plot
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
        if isa(f{j},'TIMEMATRIX') && ~israndom(f{j})
            qtmax{j} = max(max(f{j}));
            qtmin{j} = min(min(f{j}));
        elseif isa(f{j},'TIMERADIALMATRIX')
            qtmax{j}=-Inf;
            qtmin{j}=Inf;
            for kk=1:max(1,floor(length(rep)/5)):length(rep)
                qtmax{j} = max(qtmax{j},max(max(getmatrixatstep(f{j},rep(kk)))));
                qtmin{j} = min(qtmin{j},min(min(getmatrixatstep(f{j},rep(kk)))));
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
            plot(getmatrixatstep(f{j},rep(i)),varargin{1}{j},varargin{2:end});
        elseif ~ischarin('sigma',varargin) && ~ischarin('epsilon',varargin) && ~ischarin('energyint',varargin)
            plot(getmatrixatstep(f{j},rep(i)),varargin{1}{j},'FaceColor','none','EdgeColor','k',varargin{2:end});
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
