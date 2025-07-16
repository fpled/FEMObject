function mov=evolV(T,f,scanm,varargin)
% function mov=evol(T,f,varargin)
% T : TIMEMODEL
% f : repr�sente la fonction dont on
% varargin : arguments pour la fonction plot
% mov = ensemble des frames pour faire un film

if isa(f,'MULTIMATRIX') || isa(f,'PCRADIALMATRIX')
    m = length(f);
    f = TIMEMATRIX(f,T);
elseif isa(f,'PCRADIALMATRIX')
    m = length(f);
    f = TIMEMATRIX(getV(f),T);
elseif isa(f,'PCTIMEMATRIX') && isa(getvalue(f),'PCRADIALMATRIX')
    V = getV(getvalue(f));
    m = length(V);
    f = TIMEMATRIX(V,T,size(f));
elseif isa(f,'TIMEMATRIX') && isa(getvalue(f),'MULTIMATRIX')
    m = length(getvalue(f));
elseif isa(f,'cell') && isa(f{1},'TIMEMATRIX')
    m = length(f);
elseif ~isa(f,'cell')
    error('mauvais argument')
end

if isa(f,'cell')
    for k=1:length(f)
        if ~isa(f{k},'TIMEMATRIX') && ~isa(f{k},'TIMERADIALMATRIX')
            f{k} = TIMEMATRIX(f{k},T);
        end
    end
end




if isempty(scanm) || nargin<=2
    scanm = 1:m;
end

plotiter = getevolparam(T,'plotiter');
plottime = getevolparam(T,'plottime');
if plotiter || plottime
    addedplot = 1;
else
    addedplot = 0;
end

m = length(scanm);

if ischarin('subplot',varargin)
    addedplot=0;
    tempo = getcharin('subplot',varargin);
    nl=tempo(1);nc=tempo(2);
else
    nl=floor(sqrt(m+addedplot));nc=ceil((m+addedplot)/nl);
end

[t,rep]=gettevol(T);

if isa(f,'TIMEMATRIX')
    V = getvalue(f);
else
    V = f ;
end
for j=1:m
    if isa(f,'cell')
        cax{j}= [min(min(getvalue(V{j}))),max(max(getvalue(V{j})))];
    else
        cax{j}= [min(min(V{j})),max(max(V{j}))];
    end
    zli{j}= cax{j};
end

plotstep = getevolparam(T,'plotstep');
compt = 0;
for i=1:plotstep:length(t)
    compt=compt+1;
    
    execute(T.evolparam,'before');
    
    
    if plotiter
        fullsubplot(nl,nc,1);
        text(0.1,0.5,['iter ' num2str(i,'%d')],'FontSize',12)
        axis off
    elseif plottime
        fullsubplot(nl,nc,1);
        text(0.1,0.5,['time ' num2str(t(i),'%.2f') ' s'],'FontSize',12)
        axis off
    end
    
    if isa(f,'TIMEMATRIX')
        fi = f{rep(i)};
    else
        fi = f;
    end
    
    for j=1:m
        T=setevolparam(T,'caxis',cax{j});
        if ischarin('surface',varargin)
            T=setevolparam(T,'zlim',zli{j});
        end
        if ischarin('fact',varargin)
            fullsubplot(nl,nc,j+addedplot,getcharin('fact',varargin));
        else
            fullsubplot(nl,nc,j+addedplot);
        end
        
        V = fi{scanm(j)};
        if isa(V,'TIMEMATRIX')
            V = V{rep(i)};
        end
        
        manutext = getcharin('manutext',varargin);
        
        plot(V,varargin{:});
        execute(T.evolparam,'after')
        
        if ~isempty(manutext)
            modename=manutext{2};
            pos = manutext{1};
            name = manutext{2};
            
            arg = manutext(3:end);
            
            text(pos{:},[ modename '_{' num2str(scanm(j)) '}'],arg{:});
        end
        
        
    end
    
    if nargout>=1
        mov(compt) = getframe(gcf);
    end
    
end

