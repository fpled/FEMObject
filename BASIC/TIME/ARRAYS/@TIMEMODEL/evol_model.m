function mov = evol_model(T,S,varargin)
% function mov = evol_model(T,S,varargin)
% T : TIMEMODEL
% S : cell of MODEL
% varargin : arguments for function plot
% mov = movie

interpreter = getcharin('Interpreter',varargin,'latex');
leg = getcharin('legend',varargin,false);

[t,rep]=gettevol(T);

fontsize = getevolparam(T,'fontsize');
plotiter = getevolparam(T,'plotiter');
plottime = getevolparam(T,'plottime');
plotstep = getevolparam(T,'plotstep');
compt = 0;
for i=1:plotstep:length(t)
    compt=compt+1;
    execute(T.evolparam,'before');
    
    h = plot(S{i},varargin{:});
    hg = hggroup;
    set(h(:),'Parent',hg);
    if leg
        legend(hg,'$\Omega$','Interpreter',interpreter)
    end
    set(gca,'FontSize',fontsize)
    
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

return
